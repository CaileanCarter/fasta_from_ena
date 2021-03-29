#! /usr/bin/python3

import argparse
import gzip
import shutil
import urllib.request as request
import xml.etree.ElementTree as ET
from contextlib import closing
from os import path, remove, mkdir
from urllib.error import URLError
from zlib import error as zlibError

import numpy as np
import pandas as pd

f"""
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
        except AttributeError:
            unique_id = alias

        taxon = node.find('.//TAXON/SCIENTIFIC_NAME').text
        try:
            strain = node.find('.//TAXON/STRAIN').text
        except AttributeError:
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
    # df = df.loc[df["taxon"].str.startswith("Escherichia coli")]
    df = df.astype({'contigs' : int, 'length': int})

    return df


def get_FASTA(url=None, fp=None):
    with closing(request.urlopen(url)) as r, open(fp, 'wb') as f:
        shutil.copyfileobj(r, f)


def unzip_gz(zin=None, zout=None):
    try:
        with gzip.open(zin, 'rb') as f_in, open(zout, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    except zlibError:
        print(f"File could not be unzipped: {zin}")
        return
    remove(zin)


def to_excel(df=None, fp=None):
    fp = path.join(fp, "summary.xlsx")
    df.to_excel(fp)


def main(xml_fp=None, output=None):


    if not path.isdir("FASTA"):
        mkdir("FASTA")

    df = open_xml(fp=xml_fp)

    to_excel(df=df, fp=output)

    for row in df.itertuples():
        unique_id = row[0]
        url = row[6]

        unzipped = path.join(output, f"{unique_id}.fasta")
        zipped = unzipped + ".gz"
        
        try:
            get_FASTA(url=url, fp=zipped)
        except (AttributeError, URLError):
            print(f"FASTA not available for {unique_id}")
            continue
        
        unzip_gz(zin=zipped, zout=unzipped)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-i', '--input', type=str)
    parser.add_argument('-o', '--output', type=str)

    args = parser.parse_args()
    main(xml_fp=args.input, output=args.output)


if __name__ == "__main__":
    parse_arguments()
