import requests
import pandas as pd
from io import StringIO

def get_cid_from_name(drug_name):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/property/Title/CSV'
    response = requests.get(url)
    
    if response.status_code == 200:
        try:
            data_frame = pd.read_csv(StringIO(response.text))
            if not data_frame.empty:
                cid = data_frame.iloc[0]["CID"]
                return cid
        except Exception as e:
            print(f"Error processing data for drug name {drug_name}: {e}")
            return None
    else:
        print(f"Error: Unable to fetch CID for {drug_name}. Status code: {response.status_code}")
        return None

def get_all_properties_to_excel(cids, output_filename="compounds_properties.xlsx"):
    properties = [
        "Title", #Name of the molecule
        "MolecularWeight", #Weight of the molecule
        "XLogP",  # Octanol-water partition coefficient
        "ExactMass",  # Precise molecular mass
        "MonoisotopicMass",  # Mass of the most abundant isotope
        "TPSA",  # Topological Polar Surface Area
        "HBondDonorCount",  # Hydrogen bond donors
        "HBondAcceptorCount",  # Hydrogen bond acceptors
        "RotatableBondCount",  # Bond flexibility
        "HeavyAtomCount",  # Molecular backbone size
        "AtomStereoCount",  # Stereochemistry (atoms)
        "DefinedAtomStereoCount",  # Defined stereochemistry (atoms)
        "UndefinedAtomStereoCount",  # Undefined stereochemistry (atoms)
        "BondStereoCount",  # Stereochemistry (bonds)
        "DefinedBondStereoCount",  # Defined stereochemistry (bonds)
        "UndefinedBondStereoCount",  # Undefined stereochemistry (bonds)
        "CovalentUnitCount",  # Number of covalent units
        "Volume3D",  # 3D spatial extent of the molecule
        "XStericQuadrupole3D",  # 3D shape/electron density (X component)
        "YStericQuadrupole3D",  # 3D shape/electron density (Y component)
        "ZStericQuadrupole3D",  # 3D shape/electron density (Z component)
        "FeatureAcceptorCount3D",  # 3D feature acceptors
        "FeatureDonorCount3D",  # 3D feature donors
        "FeatureAnionCount3D",  # 3D feature anions
        "FeatureCationCount3D",  # 3D feature cations
        "FeatureRingCount3D",  # 3D feature rings
        "FeatureHydrophobeCount3D",  # 3D feature hydrophobes
        "ConformerModelRMSD3D",  # Variation between conformers
        "EffectiveRotorCount3D",  # Molecular flexibility
        "Fingerprint2D"  # Encoded molecular substructures
    ]


    all_data = []

    for cid in cids:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{','.join(properties)}/CSV"
        print(url)
        response = requests.get(url)

        if response.status_code == 200:
            csv_data = response.text
            data_frame = pd.read_csv(StringIO(csv_data))
            

            if "Disease" not in data_frame.columns:
                data_frame.insert(0, "Disease", "Huntington")# Replace with your disease
            else:
                data_frame["Disease"] = "Huntington"#Replace with your disease
            
            all_data.append(data_frame.iloc[0])  
            print(f"Data fetched for CID {cid}")
        else:
            print(f"Error: Unable to fetch data for CID {cid}. Status code: {response.status_code}")

    if all_data:
        combined_data_frame = pd.DataFrame(all_data)
        combined_data_frame.to_excel(output_filename, index=False)
        print(f"All data saved to {output_filename}")
    else:
        print("No data to save.")

def main(drug_names, output_filename="compounds_properties.xlsx"):
    cids = []
    for drug_name in drug_names:
        cid = get_cid_from_name(drug_name)
        if cid:
            cids.append(cid)
        else:
            print(f"Skipping {drug_name} due to missing CID.")
    
    if cids:
        get_all_properties_to_excel(cids, output_filename)
    else:
        print("No valid CIDs found. Exiting.")

#Replace with your drug names
drug_names = [
    "Tetrabenazine",
    "Tetrahydropalmatine",
    "Deutetrabenazine",
    "Valbenazine",
    "Dihydrotetrabenazine",
    "Salsolidine",
    "Carmegliptin",
    "2H-Benzo(a)quinolizin-2-ol, 2-ethyl-1,3,4,6,7,11b-hexahydro-9,10-dimethoxy-3-(2-methylpropyl)-",
    "Phellodendrine",
    "Carnegine",
    "Trifluoperazine",
    "Fluphenazine",
    "Triflupromazine",
    "Periciazine",
    "Trifluoperazine Hydrochloride",
    "Fluphenazine Hydrochloride",
    "Triflupromazine Hydrochloride",
    "Fluphenazine-N-mustard",
    "Fluphenazine Maleate",
    "Trifluoperazine dimaleate",
    "Haloperidol",
    "Haloperidol Decanoate",
    "Trifluperidol",
    "Moperone",
    "Haloperidol Lactate",
    "4-(2-(1-(4-Chlorocinnamyl)piperazin-4-yl)ethyl)benzoic acid",
    "Flazalone",
    "Moperone Hydrochloride",
    "Emetine",
    "Sulpiride",
    "Levosulpiride",
    "Sultopride",
    "Tiapride",
    "Veralipride",
    "Tiapride Hydrochloride",
    "Sultopride Hydrochloride",
    "Iodosulpride",
    "Moexipril Hydrochloride",
    "Moexipril",
    "Valbenazine tosylate",
    "N,3-dimethyl-N-[(2R)-4-(4-methylpiperidin-1-yl)butan-2-yl]benzenesulfonamide",
    "Branaplam",
    "Branaplam Hydrochloride",
    "Troriluzole",
    "Troriluzole Hydrochloride",
    "Olanzapine",
    "Quetiapine",
    "Nordazepam",
    "Clonazepam",
    "Nitrazepam",
    "Delorazepam",
    "Nimetazepam",
    "Memantine",
    "Amantadine",
    "Amantadine Hydrochloride",
    "Memantine Hydrochloride",
    "Lonafarnib",
    "Sch 59228",
    "Tricaprin",
    "Glyceryl Monostearate"

]

main(drug_names)

