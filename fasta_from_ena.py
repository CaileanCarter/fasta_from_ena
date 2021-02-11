#! /usr/bin/python3

import argparse
import gzip
import shutil
import urllib.request as request
import xml.etree.ElementTree as ET
from contextlib import closing
from os import path, remove

import numpy as np
import pandas as pd

__doc__=f"""
Download FASTA files from ENA using an XML report of a search.
Simply provide the XML file as input and specify an output directory to download FASTA files to.
An Excel file is created as a record of all the downloaded files.

Usage:
{path.basename(__file__)}.py -i [input] -o [output]
"""


def open_xml(fp=None):

    tree = ET.parse(fp)
    root = tree.getroot()

    master = {}
    for node in root:
        accession = node.attrib.get("accession")
        alias = node.attrib.get("alias")
        center_name = node.attrib.get("center_name")

        try:
            unique_id = node.find('.//WGS_SET/PREFIX').text + '0' + node.find('.//WGS_SET/VERSION').text
        except:
            unique_id = alias

        taxon = node.find('.//TAXON/SCIENTIFIC_NAME').text
        try:
            strain = node.find('.//TAXON/STRAIN').text
        except:
            strain = np.NaN

        FASTA = np.NaN
        for child in node.iterfind(".//ASSEMBLY_LINKS/ASSEMBLY_LINK/URL_LINK"):
            if child.find(".//LABEL").text == "WGS_SET_FASTA":
                FASTA = child.find(".//URL").text 

        length = np.NaN
        contigs = np.NaN

        for child in node.findall(".//ASSEMBLY_ATTRIBUTES/ASSEMBLY_ATTRIBUTE"):
            if child.find(".//TAG").text == "total-length":
                length = child.find(".//VALUE").text
            elif child.find(".//TAG").text == "count-contig":
                contigs = child.find(".//VALUE").text     
        
        master.update({
            unique_id : {
                "accession" : accession,
                "alias" : alias,
                "center_name" : center_name,
                "taxon" : taxon,
                "strain" : strain,
                "FASTA" : FASTA,
                "length" : length,
                "contigs" : contigs
            }
        })

    df = pd.DataFrame.from_dict(master, orient="index")
    df = df.loc[df["taxon"].str.startswith("Escherichia coli")]
    df = df.astype({'contigs' : int, 'length': int})

    return df


def get_FASTA(taxon=None, url=None, fp=None):
    with closing(request.urlopen(url)) as r, open(fp, 'wb') as f:
        shutil.copyfileobj(r, f)


def unzip_gz(taxon=None, zin=None, zout=None):
    with gzip.open(zin, 'rb') as f_in, open(zout, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    remove(zin)


def to_excel(df=None, fp=None):
    df.to_excel(fp)


def main(*args, xml_fp=None, output=None, **kwargs):

    df = open_xml(fp=xml_fp)

    to_excel(df=df, fp=output)

    for row in df.itertuples():
        unique_id = row[0]
        url = row[6]

        unzipped = path.join(output, f"{unique_id}.fasta")
        zipped = unzipped + ".gz"
        
        try:
            get_FASTA(taxon=unique_id, url=url, fp=zipped)
        except:
            print(f"FASTA not available for {unique_id}")
            continue
        
        unzip_gz(taxon=unique_id, zin=zipped, zout=unzipped)


parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('-i', '--input', type=str)
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()


if __name__ == "__main__":
    main(xml_fp=args.input, output=args.output)
