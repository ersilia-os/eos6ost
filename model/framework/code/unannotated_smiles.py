import sys
from typing import Callable

from utils import (
    SamplingResult,
    get_comb_of_mol_with_attachment_points,
    get_target_num_mols_for_given_mol,
    get_smiles,
)


class UnannotatedSmiles:
    def __init__(
        self, input_smile: "list[str]", total_molecules_to_generate: int, is_debug: bool
    ) -> None:
        self.input_smiles = input_smile
        self.is_debug = is_debug
        self.result: "dict[str, list[str]]" = {}
        # per input: the number of molecules to draw -> the scaffolds drawn with that number
        self.scaffolds_with_batch_size: "dict[str, dict[int, list[str]]]" = {}

        for smile in self.input_smiles:
            try:
                combinations = get_comb_of_mol_with_attachment_points((smile))
            except Exception as e:
                # e.g. a SMILES RDKit cannot read: this input gets an empty row, the others go on
                print(
                    f"[ERROR] no scaffold for {smile!r}: {type(e).__name__}: {e}",
                    file=sys.stderr,
                )
                continue
            total_combinations = sum([len(comb) for comb in combinations])

            attachment_points = 1
            molecules_with_batch_size: "dict[str, list[str]]" = {}
            for molecules in combinations:
                for mol in molecules:
                    mol_smiles = get_smiles(mol)

                    number_of_molecules_to_generate = get_target_num_mols_for_given_mol(
                        attachment_points,
                        total_combinations,
                        total_molecules_to_generate,
                    )

                    if number_of_molecules_to_generate in molecules_with_batch_size:
                        molecules_with_batch_size[
                            number_of_molecules_to_generate
                        ].append(mol_smiles)
                    else:
                        molecules_with_batch_size[number_of_molecules_to_generate] = [
                            mol_smiles
                        ]

                attachment_points += 1

            self.scaffolds_with_batch_size[smile] = molecules_with_batch_size

    def get_result(self) -> SamplingResult:
        inputs = []
        outputs = []

        for mol, res in self.result.items():
            for item in res:
                inputs.append(mol)
                outputs.append(item)

        return SamplingResult(inputs, outputs)

    def generate(
        self, sample: Callable[[list[str], int], SamplingResult]
    ) -> SamplingResult:
        # Every input is sampled on its own. When the scaffolds of all inputs shared a call, the
        # outputs were credited to the first input whose scaffold matched, so two inputs of a file
        # with the same scaffold got all the outputs for the first and an empty row for the second.
        for root, groups in self.scaffolds_with_batch_size.items():
            for batch_size, mols in groups.items():
                sampled = sample(mols, int(batch_size))
                self.result.setdefault(root, []).extend(sampled.output)

        return self.get_result()
