import os
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
            unique_sequences=self,
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

    def sample(
        self, input_smiles: "list[str]", batch_size: "int|None" = None
    ) -> SamplingResult:
        if len(input_smiles) == 0:
            return SamplingResult([], [])

        with torch.no_grad():
            sampler = self.get_sampler(batch_size)
            sampled = sampler.sample(input_smiles)

        sampled = filter_valid(sampled)
        sampled = filter_out_duplicate_molecules(sampled, is_debug=self.is_debug)

        return SamplingResult(sampled.items1, sampled.smilies)

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
