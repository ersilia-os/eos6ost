# imports
import sys
import os
import csv
import json

import click
from rdkit import Chem
from reinvent.config_parse import read_smiles_csv_file

from libinvent_sampler import LibinventSampler


def strip_atom_maps(smi):
    # REINVENT4's bond-making step leaves attachment-point atom-map numbers
    # (e.g. "[CH3:0]") on every generated SMILES; confirmed these are present in
    # 100% of outputs. They're not canonical and break naive exact-match
    # deduplication/registration downstream, so clear them before writing.
    if not smi:
        return smi
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return smi
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)


def drop_input(row, input_smi):
    # The input molecule must never be returned as one of its own "generated" outputs
    # (LibInvent regrows the original molecule from its own scaffold for some inputs:
    # 4/100 compounds in the 100-compound benchmark). Compared without stereochemistry,
    # so a stereo-stripped copy of the input is caught too. An annotated scaffold input
    # (contains "*") can never equal a complete output, so it is left alone.
    input_mol = Chem.MolFromSmiles(input_smi)
    if input_mol is None:
        return row
    input_flat = Chem.MolToSmiles(input_mol, isomericSmiles=False)
    kept = []
    for s in row:
        if s:
            mol = Chem.MolFromSmiles(s)
            if mol is not None and Chem.MolToSmiles(mol, isomericSmiles=False) == input_flat:
                continue
        kept.append(s)
    return kept

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# This arguments is reserved for testing or
# running model locally or in notebook.
is_debug = sys.argv[3] == "True" if len(sys.argv) > 3 else False

# Name of the log file.
# Only write if `is_debug` is True.
log_file = output_file + ".json"

batch_size = 1000
num_input_smiles = 0
input_smiles = None


if os.path.exists(input_file):
    input_smiles = read_smiles_csv_file(input_file, columns=0, header=True)
    num_input_smiles = len(input_smiles)

else:
    click.echo(click.style(f"[INPUT_FILE]: {input_file} doesn't exist.", fg="red"))


if not os.path.exists(os.path.dirname(os.path.abspath(output_file))):
    click.echo(
        click.style(
            f"[OUTPUT_DIR]: {os.path.dirname(output_file)} doesn't exist.", fg="red"
        )
    )

libinvent_sampler = LibinventSampler(batch_size=batch_size, is_debug=is_debug)

if is_debug:
    click.echo(click.style("Running libinvent prior", fg="white", bg="green"))

outputs, _, log_libinvent = libinvent_sampler.generate(input_smiles=input_smiles)


input_len = len(input_smiles)
output_len = len(outputs)

assert input_len == output_len

outputs = [[strip_atom_maps(s) for s in row] for row in outputs]

# drop the input molecule itself from its own row; any shortfall is padded by
# dedupe_and_pad below, like the other post-filters here (no backfill by re-generating)
outputs = [drop_input(row, smi) for row, smi in zip(outputs, input_smiles)]

# stripping atom maps can collapse two previously-distinct-looking outputs (e.g. differing
# only in which attachment point got which map number) onto the same canonical structure;
# confirmed empirically this reintroduces a small number of duplicates (~2% in one row of
# the shipped examples) that the sampler's own pre-stripping dedup couldn't have caught.
# Dedupe again here and pad any shortfall, rather than backfilling by re-generating.
def dedupe_and_pad(row, target):
    seen = set()
    deduped = []
    for s in row:
        if not s:
            continue
        if s in seen:
            continue
        seen.add(s)
        deduped.append(s)
    deduped = deduped[:target]
    return deduped + [""] * (target - len(deduped))

outputs = [dedupe_and_pad(row, batch_size) for row in outputs]

HEADER = ["smi_{0}".format(str(x).zfill(3)) for x in range(batch_size)]

with open(output_file, "w", newline="") as fp:
    csv_writer = csv.writer(fp)
    # First Row: Header
    # Second Row: Generated Smiles (Output)
    csv_writer.writerows([HEADER])
    csv_writer.writerows(outputs)


if is_debug:
    log = {
        "start": log_libinvent["start"],
        "end": log_libinvent["end"],
        "input_smiles": log_libinvent["input_smiles"],
        "total": log_libinvent["total"],
        "expected": batch_size * num_input_smiles,
    }

    with open(os.path.abspath(log_file), "w", newline="\n") as fp:
        json.dump(log, fp)
