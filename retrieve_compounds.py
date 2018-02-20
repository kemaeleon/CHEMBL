# some general packages
import json
import warnings
import logging
import pandas as pd
import sys
# loading chembl package, creating chembl python pyclient
import chembl_webresource_client
from chembl_webresource_client import *
from chembl_webresource_client.settings import Settings
from chembl_webresource_client.new_client import new_client as chembl
import xmltodict
# loading rdkit package
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import argparse
from PIL import Image, ImageFont, ImageDraw, ImageEnhance


def exists_target_for_uniprot_id(target, uniprot_id):
    try:
        if target['target_components'][0]['accession'] == uniprot_id:
            return True
    except:
        print("warning: no uniprot id was found for " + target['target_chembl_id'])
        return False
    return False


def get_target(uniprot_id):
    targets = TargetResource()
    t = targets.get(uniprot=uniprot_id)
    if 'to_bytes' in dir(t):
        sys.exit(0)
    return(t)

def has_pka_or_ic50_between(tested_compound, low_con, high_con):
    if tested_compound['standard_units'] == 'nM' and tested_compound['standard_value'] is not None \
            and float(low_con < float(tested_compound['standard_value'])) < high_con:
        return True
    return False


def has_target_confidence(tested_compound, threshold):
    chembl_assay_id = (tested_compound['assay_chembl_id'])
    details = chembl.assay.get(chembl_assay_id)
    if details['confidence_score'] >= float(threshold):
        return True
    return False


def is_assay_type(tested_compound, *types):
    if tested_compound['assay_type'] in list(types):
        return True
    return False


def depict_activity(a_s, rank_low, rank_high, out_name, is_drug=False):
    selected_ranks = []
    if rank_low < 0:
        rank_high = len(a_s)
        rank_low = rank_high + rank_low
    for ii in range(rank_low, rank_high):
        z = a_s.iloc[ii]
        base_label = str(z['standard_value']) + str(z['standard_units'])
        if not is_drug:
            label = base_label
        else:
            label = str(z['drug_name']) + " " + base_label
        selected_ranks.append((z['canonical_smiles'], label))
    molecules = []
    for a in selected_ranks:
        m = Chem.MolFromSmiles(a[0])
        AllChem.Compute2DCoords(m)
        m.SetProp('_Name', a[1])
        molecules.append(m)
    img = Draw.MolsToGridImage(molecules, molsPerRow=3, subImgSize=(200, 200),
                               legends=[m.GetProp('_Name') for m in molecules])
    cropped = img.crop((0,-34,img.width,img.height))
    draw = ImageDraw.Draw(cropped)
    draw.rectangle((0, -34, img.width, 34), fill="white")
    draw.text((100, 20), str(translate_label(out_name)), font=ImageFont.truetype("arial.ttf", 14),
            fill='#000000')
    canvas = Image.new('RGBA', cropped.size, (255,255,255,255))
    canvas.paste(cropped, mask=cropped)
    canvas.save(out_name)


def translate_label(i_out_name):
    parts = i_out_name.split('.')[0].split('_')
    if len(parts) == 4:
        return "TARGET ID: {} {} {}".format(parts[0], parts[2], parts[3])
    else:
        return "Target ID: {} ".format(parts[0])


def get_known_drugs(i_sorted_active_selection, ichembl_id):
    targets = TargetResource()
    drugs = targets.approved_drugs(ichembl_id)
    i_known_drugs_selection = []
    for ix in i_sorted_active_selection:
        for iy in drugs:
            if iy['chemblId'] == ix['molecule_chembl_id']:
                tmp = ix
                tmp['drug_name'] = iy['name']
                i_known_drugs_selection.append(tmp)
    return i_known_drugs_selection


def retrieve_act(tested_compounds, n_struct_depict, target_confidence, var_string):
    """retrieve compound structure and bioactivity"""
    active_selection = [a for a in tested_compounds if has_pka_or_ic50_between(a, 1e-3, 1e7) and
                        a['assay_type'] in ['B', 'F'] and has_target_confidence(a, target_confidence)]
    sorted_active_selection = sorted(active_selection, key=lambda k: float(k['standard_value']))
    df = pd.DataFrame(sorted_active_selection)
    uniq_active_selection = df.drop_duplicates(subset=['canonical_smiles'])
    print("...number of compounds that matched our criteria " + str(len(uniq_active_selection)))
    depict_activity(uniq_active_selection, 0, min(n_struct_depict, len(uniq_active_selection)),
                    "{}_best_actives.png".format(var_string))
    depict_activity(uniq_active_selection, -min(len(active_selection), n_struct_depict), 0,
                    "{}_least_actives.png".format(var_string))


def retrieve_drugs(sorted_active_selection, chemblId, var_string):
    """retrieve compound structure and bioactivity of known drugs"""
    known_drug_selection = get_known_drugs(sorted_active_selection,chemblId)
    if len(known_drug_selection) == 0:
        print("...no approved drugs with CHEMBL bioactivity measurements could be found")
        sys.exit(0)
    df_drugs = pd.DataFrame(known_drug_selection)
    df_drugs = df_drugs.drop_duplicates(subset=['canonical_smiles'])
    print("...number of structurally uniq drugs for target " + str(len(df_drugs)))
    depict_activity(df_drugs, 0, len(df_drugs),
                    "{}_drugs.png".format(var_string), True)


def main():
    parser = argparse.ArgumentParser(description='Retrieve structure and activity of compounds stored in CHEMBL')
    parser.add_argument('--uniprot_id', dest='uniprot_id', type=str, default='P00519',
                        help='uniprot_id of your target')
    parser.add_argument('--n_struct_pic', type=int, default=6, help='number of top active and top inactive compounds \
                        to depict, this is ignored if --drugs option is selected')
    parser.add_argument('--drugs', type=bool, dest='drugs', default=False, help='activities of known drugs only')
    parser.add_argument('--target_confidence', type=int, default=9, help='minimum accepted target confidence')
    args = parser.parse_args()
    target = get_target(args.uniprot_id)
    var_string = ''
    if not args.drugs:
        var_string = "{}_{}".format(target['chemblId'], args.n_struct_pic)
    else:
        var_string = "{}".format(target['chemblId'])
    print("...retrieved target " + target['chemblId'] + " " +
          target['description'] + " "  + target['synonyms'] )
    tested_compounds = chembl.activity.filter(target_chembl_id=target['chemblId'])
    if not args.drugs:
        retrieve_act(tested_compounds, args.n_struct_pic, args.target_confidence, var_string)
    else:
        retrieve_drugs(tested_compounds, target['chemblId'], var_string)


if __name__ == '__main__':
    main()
