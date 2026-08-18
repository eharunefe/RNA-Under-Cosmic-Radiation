import numpy as np
import pandas as pd
from scipy import stats
import urllib.request
import os
import json
import webbrowser


df = pd.read_csv("cumulative.csv")


clean_df = df.dropna(subset=["SYMBOL"]).copy()
samples = [c for c in clean_df.columns if c.startswith("GSM")]





genenames = clean_df.groupby("SYMBOL")["GENENAME"].first()
samples_df = clean_df.groupby("SYMBOL")[samples].mean()



samples_df.insert(0, "GENE NAME", genenames)


print("Tablo boyutu")
print("genler " + str(samples_df.shape[0]))
print("nümuneler " + str(samples_df.shape[1] - 1))


print(samples_df.head())


cont = ["GSM506087", "GSM506088", "GSM506089", "GSM506090", "GSM506091"]
rady = ["GSM506092", "GSM506093", "GSM506094", "GSM506095", "GSM506096"]
micr = ["GSM506097", "GSM506098", "GSM506099", "GSM506100", "GSM506101"]
deep = ["GSM506102", "GSM506103", "GSM506104", "GSM506115", "GSM506116"]


samples_df["Earth Normal"] = samples_df[cont].mean(axis=1)
samples_df["Space w Radiation"] = samples_df[deep].mean(axis=1)



samples_df["Gap(c-d)"] = (

    samples_df["Space w Radiation"] - samples_df["Earth Normal"]

)



t_stat, p_val = stats.ttest_ind(

    samples_df[deep], samples_df[cont], axis=1, equal_var=False

)


samples_df["p value"] = p_val
samples_df["pValue w -Log10"] = -np.log10(samples_df["p value"])

real_samples_df = samples_df[samples_df["p value"] < 0.05]
real_samples_df = real_samples_df[samples_df["Gap(c-d)"] > 1]

tops = real_samples_df.sort_values(by="Gap(c-d)", ascending=False)[

    ["GENE NAME", "Gap(c-d)", "p value", "pValue w -Log10"]

].head(5)

genes = real_samples_df.sort_values(by="Gap(c-d)", ascending=False)[

    ["GENE NAME", "Gap(c-d)", "p value", "pValue w -Log10"]

].head(samples_df.shape[0])

print(genes)

genes = genes.index.tolist()

protein = {}
def sentez(gen):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=(gene_exact:{gen})+AND+(organism_id:9606)+AND+(reviewed:true)&format=fasta&size=1"
    ist = urllib.request.Request(
        url, headers={"User-Agent": "SpaceBio/Python"}
    )


    try:
        with urllib.request.urlopen(ist) as turn:
            context = turn.read().decode("utf-8")

        rw = context.strip().split("\n")
        uniprotID = rw[0].split("|")[1]
        aminos = "".join(rw[1:])

        return uniprotID, aminos
    except:
        return None, None



for gen in genes:
    uid, dizi = sentez(gen)
    if uid != None:
        protein[gen] = {"id": uid, "dizi": dizi, "lenght": len(dizi)}
    else:
        protein[gen] = {"id": uid, "dizi": dizi, "lenght": 0}
    print(f"{list(protein.keys()).index(gen) + 1}/{len(genes)}")


weak = {}
for gen in protein:
    weak[gen] = 0
    try:
        for c in protein[gen]["dizi"]:
            if c == "C":
                weak[gen] += 2
            elif c == "M":
                weak[gen] += 1
            elif c in ["W", "Y", "H"]:
                weak[gen] += 0.5

        weak[gen] = weak[gen] / len(protein[gen]["dizi"])
    except:
        pass


    real_samples_df.loc[gen, "Weakness"] = weak[gen] * 100
    real_samples_df.loc[gen, "Paradox"] = weak[gen] * 2**real_samples_df.loc[gen, "Gap(c-d)"]

net_df = real_samples_df.sort_values(by="Paradox", ascending=False)[

    ["GENE NAME", "Gap(c-d)", "Weakness", "Paradox"]

]
net_df = net_df[net_df["Paradox"] > 0]
top_gens = net_df.index.tolist()
print(net_df)


pdb = {}
for top_gen in top_gens:
    top_uid = protein[top_gen]["id"]

    api_ist = urllib.request.Request(
        f"https://alphafold.ebi.ac.uk/api/prediction/{top_uid}",
        headers={"User-Agent": "SpaceBio/Python"},
    )
    with urllib.request.urlopen(api_ist) as cvp:
        pdb_url = json.loads(cvp.read().decode("utf-8"))[0]["pdbUrl"]

    pdb_ist = urllib.request.Request(
        pdb_url, headers={"User-Agent": "SpaceBio/Python"}
    )
    with urllib.request.urlopen(pdb_ist) as cvp:
        pdb[top_gen] = cvp.read().decode("utf-8")






# 1. Derin Uzay (Deep Space) Verisinden Çıkan 6 Gen İçin Açıklamalar
sums = {
    "CCL2": (
        "Recruits cleanup immune cells to damage sites, but radiation cuts the"
        " fragile bridges holding it together."
    ),
    "TMEM176A": (
        "Controls cellular ion balance and immune tolerance, but its membrane"
        " structure degrades under cosmic radiation."
    ),
    "HDC": (
        "Produces histamine to trigger emergency inflammation, but radiation"
        " breaks its sensitive catalytic core."
    ),
    "PPA1": (
        "Provides essential energy for DNA repair and metabolism, but loses"
        " function when oxidized by free radicals."
    ),
    "ADAMDEC1": (
        "Helps remodel damaged tissue and defense signals, but is packed with"
        " sulfur bonds easily destroyed by radiation."
    ),
    "UTS2": (
        "Regulates cardiovascular stress, but its cyclic shape is instantly"
        " unlocked and inactivated by radiation."
    ),
}

# 1. Panelleri Oluştur (Sembol, Tam Gen Adı, Paradoks Skoru ve Açıklama ile)
polidivs = "\n    ".join(
    [
        f"""<div class="panel-wrap">
        <div class="info-box">
            <div class="gene-title">
                {g} 
                <span class="full-name">({net_df.loc[g, 'GENE NAME']})</span>
                <span class="score">Paradox: {round(net_df.loc[g, 'Paradox'], 2)}</span>
            </div>
            <div class="gene-desc">{sums.get(g, 'Heavily produced during spaceflight stress, but highly vulnerable to radiation damage.')}</div>
        </div>
        <div id="view_{g}" class="panel"></div>
    </div>"""
        for g in top_gens
    ]
)

js_cagrilari = "\n        ".join(
    [f'ciz("view_{g}", {json.dumps(pdb.get(g, ""))});' for g in top_gens]
)

# 2. 3x2 Grid HTML Şablonu
html = f"""<!DOCTYPE html>
<html>
<head>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{
            margin: 0;
            background: black;
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            grid-template-rows: repeat(2, 1fr);
            width: 100vw;
            height: 100vh;
            overflow: hidden;
            box-sizing: border-box;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        }}
        .panel-wrap {{
            position: relative;
            width: 100%;
            height: 100%;
            border: 1px solid #1a1a1a;
            box-sizing: border-box;
        }}
        .panel {{
            width: 100%;
            height: 100%;
        }}
        .info-box {{
            position: absolute;
            top: 8px;
            left: 8px;
            right: 8px;
            background: rgba(10, 10, 10, 0.85);
            backdrop-filter: blur(4px);
            border: 1px solid #27272a;
            border-radius: 5px;
            padding: 6px 10px;
            z-index: 10;
            pointer-events: none;
        }}
        .gene-title {{
            color: #ffffff;
            font-size: 13px;
            font-weight: 700;
            margin-bottom: 2px;
            display: flex;
            align-items: center;
            flex-wrap: wrap;
            gap: 4px;
        }}
        .full-name {{
            color: #94a3b8;
            font-size: 11px;
            font-weight: 400;
            font-style: italic;
        }}
        .score {{
            color: #f59e0b;
            font-size: 11px;
            font-weight: 600;
            margin-left: auto;
        }}
        .gene-desc {{
            color: #d4d4d8;
            font-size: 10.5px;
            line-height: 1.35;
            margin-top: 2px;
        }}
    </style>
</head>
<body>
    {polidivs}

    <script>
        function ciz(id, pdbVerisi) {{
            if (!pdbVerisi) return;
            let el = document.getElementById(id);
            let v = $3Dmol.createViewer(el, {{backgroundColor: "black"}});
            v.addModel(pdbVerisi, "pdb");
            v.setStyle({{}}, {{cartoon: {{color: 'grey'}}}});
            v.addStyle({{resn: 'CYS'}}, {{stick: {{color: '#FF0400'}}, sphere: {{scale: 0.35, color: '#FF0400'}}}});
            v.addStyle({{resn: 'MET'}}, {{stick: {{color: '#FF9500'}}, sphere: {{scale: 0.35, color: '#FF9500'}}}});
            v.addStyle({{resn: ['TRP', 'TYR', 'HIS']}}, {{stick: {{color: '#FFFB00'}}}});
            v.zoomTo();
            v.render();
        }}

        {js_cagrilari}
    </script>
</body>
</html>"""

with open("model.html", "w", encoding="utf-8") as f:
    f.write(html)

webbrowser.open("file://" + os.path.realpath("model.html"))