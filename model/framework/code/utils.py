from itertools import combinations
from math import ceil

import click
from reinvent.models.model_factory.sample_batch import SampleBatch
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold


class SamplingResult:
    """A class to save the sampling result."""

    def __init__(self, input_items: "list[str]", output_items: "list[str]") -> None:
        self.input = input_items
        self.output = output_items

    def append(self, other: "SamplingResult"):
        self.input.extend(other.input)
        self.output.extend(other.output)
        pass


def are_smiles_same(smiles1: str, smiles2: str) -> bool:
    """To check whether given smiles are same."""
    # Parse SMILES strings to obtain RDKit molecule objects
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    # Check if both SMILES strings are valid
    if mol1 is None or mol2 is None:
        return False  # Unable to parse one or both SMILES strings

    # Generate Morgan fingerprints for each molecule
    fp1 = Chem.AllChem.GetMorganFingerprint(mol1, 2)
    fp2 = AllChem.GetMorganFingerprint(mol2, 2)

    # Compare fingerprints to assess similarity
    try:
        similarity = float(DataStructs.TanimotoSimilarity(fp1, fp2))
    except ValueError:
        similarity = -1
    # Determine whether the molecules are considered the same
    return similarity == 1.0


def filter_out_duplicate_molecules(sampled: SampleBatch, is_debug: bool) -> SampleBatch:
    """Filter out duplicate molecules from the sampled molecules.
    It also remove the output molecules if it is similar to input molecules.

    `sampled.items1` contains input smiles.
    `sampled.smilies` contains output smiles.
    """

    seen = {}
    items1 = []
    items2 = []
    states = []
    smilies = []
    nll_indices = []

    for item1, item2, smile, nll_index, state in zip(
        sampled.items1,
        sampled.items2,
        sampled.smilies,
        range(len(sampled.items1)),
        sampled.states,
    ):
        seen[item1] = item1

        if smile in seen:
            if is_debug:
                click.echo(
                    click.style(
                        f"Removing {smile}, as it is a duplicate entry.", fg="yellow"
                    )
                )
            continue

        seen[smile] = smile

        items1.append(item1)
        items2.append(item2)
        smilies.append(smile)
        states.append(state)
        nll_indices.append(nll_index)

    return SampleBatch(
        items1=items1,
        items2=items2,
        states=states,
        smilies=smilies,
        nlls=sampled.nlls[nll_indices],
    )


def pad_smiles(
    sampled: SamplingResult, input_smiles: "list[str]", target_length=100
) -> "list[str]":
    """For a given input smiles, it is not possible to get the expected
    number (target_length) of output smiles (samples). This will cause
    problems in the downstream process. To mitigate this, the function
    pads the output with empty strings to match the target_length.
    """

    output_smiles: dict[str, list[str]] = {}

    # Initialize output_smiles with empty lists for each input smile
    output_smiles = {smile: [] for smile in input_smiles}

    # Smiles with similar tanimoto.
    similar_smiles = {}

    # Populate output_smiles with sampled smiles
    for idx, seq in enumerate(sampled.input):
        if seq in output_smiles:
            output_smiles[seq].append(sampled.output[idx])
        else:
            # REINVENT4 do some modification in the input smiles. Like if the
            # input smile is `CC(=O)Oc1ccccc1C(O)=O` then it will convert it
            # to `CC(=O)Oc1ccccc1C(=O)O`. Both of them are same, however both
            # strings are not same. This else condition will take care of this
            # edge case.
            if seq in similar_smiles:
                output_smiles[similar_smiles[seq]].append(sampled.output[idx])
            else:
                # Try to find similar smiles in input_smiles.
                for smile in input_smiles:
                    if are_smiles_same(smile, seq):
                        output_smiles[smile].append(sampled.output[idx])
                        similar_smiles[seq] = smile
                        break

    output = []

    # Construct the output list
    for smile in input_smiles:
        output_smile = output_smiles[smile][:target_length]
        output_smile_length = len(output_smile)
        if output_smile_length == target_length:
            output.extend(output_smile)
        else:
            output.extend(output_smile)
            padding = [""] * (target_length - output_smile_length)
            output.extend(padding)

    return output


def make_list_into_lists_of_n(lst: "list[str]", n: int) -> "list[list[str]]":
    """This function splits a list into n parts of equal size."""

    length = len(lst)

    if length % n != 0:
        # If length is not divisible by n, then
        # We don't have expected output.
        # Hence throwing error.
        raise Exception(f"{length} is not divisible by {n}")

    part_length = length // n

    output = []

    for start in range(0, length, part_length):
        end = start + part_length
        output.append(lst[start:end])

    return output


def get_smiles(mol):
    return Chem.MolToSmiles(mol)


def get_mol(smiles):
    return Chem.MolFromSmiles(smiles)


def get_scaffold(mol):
    return MurckoScaffold.GetScaffoldForMol(mol)


def get_idxs_of_carbon_for_new_bond(mol):
    """Return indices of all carbon atoms available for new bond.

    Technically, it returns the indices of Carbon atoms that have
    at least one bond with hydrogen atom.
    """

    carbon_indices = []

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            bonds = atom.GetBonds()
            num_bonds = sum([b.GetBondTypeAsDouble() for b in bonds])
            if num_bonds < 4:
                carbon_indices.append(atom.GetIdx())

    return carbon_indices


def get_scaffold_and_attachment_points(mol):
    scaffold = get_scaffold(mol)

    if Chem.MolToSmiles(scaffold) == "":
        # If we have a valid scaffold then we are using it.
        scaffold = Chem.Mol(mol)

    return (scaffold, get_idxs_of_carbon_for_new_bond(scaffold))


def get_mol_after_adding_attachment_points_at(mol, at):
    connecting_atom = Chem.Atom("*")
    mutable_copy = Chem.RWMol(mol)

    for attachment_idx in at:
        connection_idx = mutable_copy.AddAtom(connecting_atom)
        mutable_copy.AddBond(attachment_idx, connection_idx, Chem.BondType.SINGLE)

    _mol = mutable_copy.GetMol()
    AllChem.Compute2DCoords(_mol)

    return _mol


def attach_num_to_attachment_points(mol):
    mol_str = get_smiles(mol)
    smiles = ""
    count = 1
    for char in mol_str:
        if char == "*":
            smiles += f"[*:{count}]"
            count += 1
        else:
            smiles += char

    return get_mol(smiles)


def get_target_num_mols_for_given_mol(
    attachment_points_using, total_combinations, target=1000
):
    target_for_current_mol = target // total_combinations

    muliplier = 1.5
    if attachment_points_using == 2:
        muliplier = 2
    if attachment_points_using == 3:
        muliplier = 2.5

    # Usually there is a lot of loss when we try to generate new molecules using libinvent
    # Like, if we ask to generate 100 new molecules libinvent could generate around 60 or 70.
    # And if the input molecule has more attachment point then it could generate around 40 or 30.
    # Therefore we are asking to generate more in order to minimise the loss.
    return ceil(target_for_current_mol * muliplier)


def filter_duplicate_molecules(molecules):
    _mols = [get_smiles(mol) for mol in molecules]
    seen = {}
    filtered_list = []

    for m in _mols:
        if m in seen:
            continue
        seen[m] = m
        filtered_list.append(get_mol(m))

    return filtered_list


def get_comb_of_mol_with_attachment_points(smile):
    mol = get_mol(smile)
    scaffold, at = get_scaffold_and_attachment_points(mol)

    if len(at) > 5:
        # Limiting max attachment point to 5
        at = at[:5]

    combs = [combinations(at, 1), combinations(at, 2), combinations(at, 3)]

    # At index 0 we have all molecules with one attachment points.
    # At index 1 we have all molecules with two attachment points.
    # At index 2 we have all molecules with three attachment points.
    molecules_list = []
    attachment_point = 0

    for comb in combs:
        molecules_list.append([])

        for current_at in comb:
            s = Chem.Mol(scaffold)
            m = get_mol_after_adding_attachment_points_at(s, current_at)
            # m = attach_num_to_attachment_points(m)
            molecules_list[attachment_point].append(m)

        attachment_point += 1

    return [filter_duplicate_molecules(molecules) for molecules in molecules_list]


def print_details(combinations):
    total = sum([len(a) for a in combinations])

    attachment_points = 1
    out = {}

    for molecules in combinations:
        for mol in molecules:
            s = get_smiles(mol)
            t = get_target_num_mols_for_given_mol(attachment_points, total)

            if t in out:
                out[t].append(s)
            else:
                out[t] = [s]
        attachment_points += 1
    print(out)
