from typing import Callable, TypedDict

from utils import (
    SamplingResult,
    get_comb_of_mol_with_attachment_points,
    get_target_num_mols_for_given_mol,
    get_smiles,
    are_smiles_same,
)


class AnnotatedSmilesDict(TypedDict):
    root: str
    mol: str


class UnannotatedSmiles:
    def __init__(
        self, input_smile: "list[str]", total_molecules_to_generate: int, is_debug: bool
    ) -> None:
        self.input_smiles = input_smile
        self.is_debug = is_debug
        self.result: "dict[str, list[str]]" = {}
        self.annotated_smiles_with_batch_size: "dict[str, list[AnnotatedSmilesDict]]" = {}
        self.annotation_to_root_map: "dict[str, str]" = {}

        for smile in self.input_smiles:
            combinations = get_comb_of_mol_with_attachment_points((smile))
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

            for batch_size, molecules in molecules_with_batch_size.items():
                entries: "list[AnnotatedSmilesDict]" = [
                    {"root": smile, "mol": mol} for mol in molecules
                ]

                if batch_size in self.annotated_smiles_with_batch_size:
                    self.annotated_smiles_with_batch_size[batch_size].extend(entries)
                else:
                    self.annotated_smiles_with_batch_size[batch_size] = entries

    def find_root(
        self, molecules: "list[AnnotatedSmilesDict]", smile: str
    ) -> str | bool:
        if smile in self.annotation_to_root_map:
            return self.annotation_to_root_map[smile]

        for mol in molecules:
            if are_smiles_same(smile, mol["mol"]):
                self.annotation_to_root_map[smile] = mol["root"]
                return self.annotation_to_root_map[smile]
        return False

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
        for batch_size, molecules in self.annotated_smiles_with_batch_size.items():
            batch = int(batch_size)
            mols = [d["mol"] for d in molecules]
            sampled = sample(mols, batch)

            for mol, res in zip(sampled.input, sampled.output):
                root = self.find_root(molecules, mol)

                if root is False:
                    if self.is_debug:
                        print(f"unable to find root for {mol}")
                    continue

                if root in self.result:
                    self.result[root].append(res)
                else:
                    self.result[root] = [res]

        return self.get_result()
