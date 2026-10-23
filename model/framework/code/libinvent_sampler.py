import os
import sys
import time

import click
import torch
from reinvent.runmodes import samplers, create_adapter
from reinvent.runmodes.dtos import ChemistryHelpers
from reinvent.chemistry import TransformationTokens, Conversions
from reinvent.chemistry.library_design import BondMaker, AttachmentPoints
from reinvent.runmodes.samplers.run_sampling import filter_valid

from unannotated_smiles import UnannotatedSmiles
from utils import (
    filter_out_duplicate_molecules,
    pad_smiles,
    make_list_into_lists_of_n,
    SamplingResult,
)


class LibinventSampler:
    ROOT = os.path.dirname(os.path.abspath(__file__))
    CHECKPOINT = os.path.join(ROOT, "..", "..", "checkpoints")
    LIBINVENT_MODEL = os.path.realpath(os.path.join(CHECKPOINT, "libinvent.prior"))

    def __init__(self, batch_size: int, is_debug=False):
        self.is_debug = is_debug
        self.batch_size = batch_size
        self.chemistry = ChemistryHelpers(
            Conversions(), BondMaker(), AttachmentPoints()
        )

        # Constants
        self.temperature = 1.0
        self.unique_sequences = True
        self.isomeric = True
        self.sample_strategy = "beamsearch"
        self.tokens = TransformationTokens()
        self.randomize_smiles = True

        # Creating adapter
        libinvent_agent, _, _ = create_adapter(
            self.LIBINVENT_MODEL, "inference", torch.device("cpu")
        )

        self.libinvent_agent = libinvent_agent

    def get_sampler(self, batch_size: "int|None" = None):
        return samplers.LibinventSampler(
            self.libinvent_agent,
            batch_size=self.batch_size if batch_size is None else batch_size,
            sample_strategy=self.sample_strategy,
            isomeric=self.isomeric,
            randomize_smiles=self.randomize_smiles,
            unique_sequences=self.unique_sequences,
            chemistry=self.chemistry,
            temperature=self.temperature,
            tokens=self.tokens,
        )

    def annotated_and_unannotated_smiles(
        self, input_smiles: "list[str]"
    ) -> "tuple[list[str], list[str]]":
        annotated = []
        unannotated = []

        for smile in input_smiles:
            if "*" in smile:
                annotated.append(smile)
            else:
                unannotated.append(smile)

        return annotated, unannotated

    def can_encode(self, sampler, scaffold: str) -> bool:
        """Whether the scaffold vocabulary of LibInvent can represent the scaffold.

        REINVENT silently drops a scaffold that has a token outside its vocabulary
        (phosphorus, selenium, a protonated nitrogen, ...) and raises an AssertionError
        when none is left, which used to take down the whole input file. The scaffold is
        tokenised the way REINVENT does it (randomised, attachment numbers removed).
        """
        vocabulary = self.libinvent_agent.get_vocabulary()
        try:
            randomized = sampler._get_randomized_smiles([scaffold])[0]
            clean = self.chemistry.attachment_points.remove_attachment_point_numbers(
                randomized
            )
            tokens = vocabulary.scaffold_tokenizer.tokenize(clean)
            return vocabulary.scaffold_vocabulary.encode(tokens) is not None
        except Exception:
            return False

    def sample_group(self, sampler, input_smiles: "list[str]") -> SamplingResult:
        with torch.no_grad():
            sampled = sampler.sample(input_smiles)

        sampled = filter_valid(sampled)
        sampled = filter_out_duplicate_molecules(sampled, is_debug=self.is_debug)

        return SamplingResult(sampled.items1, sampled.smilies)

    def sample(
        self, input_smiles: "list[str]", batch_size: "int|None" = None
    ) -> SamplingResult:
        if len(input_smiles) == 0:
            return SamplingResult([], [])

        sampler = self.get_sampler(batch_size)
        encodable = [s for s in input_smiles if self.can_encode(sampler, s)]
        if len(encodable) < len(input_smiles):
            print(
                f"[WARN] skipped {len(input_smiles) - len(encodable)} of "
                f"{len(input_smiles)} scaffolds the model cannot represent",
                file=sys.stderr,
            )
        input_smiles = encodable
        if len(input_smiles) == 0:
            return SamplingResult([], [])

        try:
            return self.sample_group(sampler, input_smiles)
        except Exception as e:
            # A scaffold the model chokes on must not blank the other scaffolds of the group:
            # redo them one at a time and skip the failing ones.
            print(
                f"[ERROR] sampling a group of {len(input_smiles)} scaffolds failed "
                f"({type(e).__name__}: {e}); retrying one at a time",
                file=sys.stderr,
            )
            result = SamplingResult([], [])
            for scaffold in input_smiles:
                try:
                    result.append(self.sample_group(sampler, [scaffold]))
                except Exception as e_single:
                    print(
                        f"[ERROR] no outputs for scaffold {scaffold!r}: "
                        f"{type(e_single).__name__}: {e_single}",
                        file=sys.stderr,
                    )
            return result

    def generate(
        self, input_smiles: "list[str]"
    ) -> "tuple[list[list[str]], list[str], dict]":
        start_time = time.time()
        num_input_smiles = len(input_smiles)

        if self.is_debug:
            click.echo(
                click.style(
                    f"Starting sampling at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}",
                    fg="green",
                )
            )
            click.echo(f"Total input smiles: {num_input_smiles}")

        annotated, unannotated = self.annotated_and_unannotated_smiles(input_smiles)
        result = self.sample(annotated)
        unannotated = UnannotatedSmiles(unannotated, self.batch_size, self.is_debug)
        result_2 = unannotated.generate(self.sample)
        result.append(result_2)

        end_time = time.time()

        if self.is_debug:
            click.echo(
                click.style(
                    f"Time taken in seconds: {int(end_time - start_time)}",
                    fg="green",
                )
            )

        flatten_outputs = pad_smiles(
            result, input_smiles=input_smiles, target_length=self.batch_size
        )

        output_smiles = make_list_into_lists_of_n(flatten_outputs, num_input_smiles)

        log = {
            "start": 0,
            "end": 0,
            "input_smiles": input_smiles,
            "total": 0,
            "expected": 0,
        }

        if self.is_debug:
            total_smiles = 0
            expected_num_smiles = num_input_smiles * self.batch_size

            for smile in flatten_outputs:
                if smile != "":
                    total_smiles += 1

            click.echo(
                click.style(
                    f"Total unique smiles generated: {total_smiles}, Expected: {expected_num_smiles}, Loss: {expected_num_smiles - total_smiles}",
                    fg="green",
                )
            )

            log = {
                "start": start_time,
                "end": end_time,
                "input_smiles": input_smiles,
                "total": total_smiles,
                "expected": expected_num_smiles,
            }

        return (output_smiles, flatten_outputs, log)
