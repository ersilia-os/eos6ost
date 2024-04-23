from utils import SamplingResult, get_all_comb_of_mol_with_attachment_points


class ProcessUnannotatedSmiles:
    def __init__(self, input_smile) -> None:
        self.input_smiles = input_smile
        self.output: "dict[str, list[str]]" = {}
        self.annotated_smiles_with_batch_size = {}

        # use get_all_comb_of_mol_with_attachment_points to generate all the
        # smiles with annotation.

    def get_annotated_smiles_for_batch_size(self, batch_size="int"):
        if batch_size in self.annotated_smiles_with_batch_size:
            return self.annotated_smiles_with_batch_size[batch_size]
        return []

    def set_intermediate_result(self, result: SamplingResult):
        for input, output in zip(result.input, result.output):
            if input in self.output:
                self.output[input].append(output)
            else:
                self.output[input] = [output]

    def get_sampled_result(self):
        output = self.output.values()
        length = len(output)
        input = [self.input_smiles for i in range(length)]

        return SamplingResult(input, output)
