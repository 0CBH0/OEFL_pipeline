import sys, re, copy, os, codecs, gzip, pysam, tables, Levenshtein
import pandas as pd
import numpy as np
import subprocess as sp
from collections import defaultdict
from scipy import sparse
from Bio import SeqIO

class as_info:
    def __init__(self):
        self.type = ""
        self.trans = []

class exon_info:
    def __init__(self):
        self.start = 0
        self.end = 0
        self.cds_s = 0
        self.cds_e = 0
        self.cds_sm = 0
        self.cds_em = 0

class transcript_info:
    def __init__(self):
        self.type = ""
        self.exon = []
        self.start = 0
        self.end = 0
        self.utr3 = 0
        self.utr5 = 0
        self.body = 0
        self.cds = ""
        self.id = ""

#pwk_genome = SeqIO.to_dict(SeqIO.parse("/home/cbh/library/pwk_filter.fa", "fasta"))
genome = SeqIO.to_dict(SeqIO.parse("/home/cbh/library/refdata-gex-mm10-2020-A/fasta/genome.fa", "fasta"))
chr_list = list(set(["chr"+str(x) for x in list(range(1, 50))+["X", "Y"]]).intersection(set(genome.keys())))
ref = pd.read_table("/home/cbh/library/refdata-gex-mm10-2020-A/genes/genes.gtf", comment="#", header=None, dtype={0:str})
ref.columns = ["seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
ref_cds = ref[(ref["type"] == "CDS") & (ref["seq_id"].isin(chr_list))].copy(deep=True)
ref_cds["trans"] = list(map(lambda x: re.findall(r'transcript_id "(.*?)";', x)[0], ref_cds["attributes"]))
ref_cds["gene"] = list(map(lambda x: re.findall(r'gene_id "(.*?)";', x)[0], ref_cds["attributes"]))
ref_cds = ref_cds[["seq_id", "start", "end", "strand", "trans", "gene"]].sort_values(by=["seq_id", "start", "end", "strand", "trans", "gene"])
#ref_cds.drop_duplicates(subset=["seq_id", "start", "end", "strand", "trans", "gene"], keep="first", inplace=True)
for i in range(ref_cds.shape[0]):
    if ref_cds.iloc[i, 2] < ref_cds.iloc[i, 1]:
        t = ref_cds.iloc[i, 2]
        ref_cds.iloc[i, 2] = ref_cds.iloc[i, 1]
        ref_cds.iloc[i, 1] = t

ref_exon = ref[(ref["type"] == "exon") & (ref["seq_id"].isin(chr_list))].copy(deep=True)
ref_exon["trans"] = list(map(lambda x: re.findall(r'transcript_id "(.*?)";', x)[0], ref_exon["attributes"]))
ref_exon["gene"] = list(map(lambda x: re.findall(r'gene_id "(.*?)";', x)[0], ref_exon["attributes"]))
ref_exon["type"] = list(map(lambda x: re.findall(r'transcript_type "(.*?)";', x)[0], ref_exon["attributes"]))
ref_exon = ref_exon[["seq_id", "start", "end", "strand", "trans", "gene", "type"]].sort_values(by=["seq_id", "start", "end", "strand", "trans", "gene"])
#ref_exon.drop_duplicates(subset=["seq_id", "start", "end", "strand", "trans", "gene"], keep="first", inplace=True)
exon_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(exon_info))))
for index, row in ref_exon.iterrows():
    chr, start, end, strand, trans, gene, type = list(row)
    exon_dict[chr][gene+":"+strand][str(start)+"-"+str(end)][trans].start = min(start, end)
    exon_dict[chr][gene+":"+strand][str(start)+"-"+str(end)][trans].end = max(start, end)

for chr, genes in exon_dict.items():
    for gene, exons in genes.items():
        strand = gene[-1]
        ref_cds_gene = ref_cds[ref_cds["gene"] == gene[:-2]]
        for k, exon in exons.items():
            for trans, e in exon.items():
                ref_cds_trans = ref_cds_gene[(ref_cds_gene["trans"] == trans) & (ref_cds_gene["start"] >= e.start) & (ref_cds_gene["end"] <= e.end)]
                if ref_cds_trans.shape[0] > 0:
                    if ref_cds_trans.shape[0] > 1:
                        print("???")
                    e.cds_sm = ref_cds_trans.iloc[0, 1]
                    e.cds_em = ref_cds_trans.iloc[0, 2]
                elif strand == "+":
                    cds = ref_cds_gene[(ref_cds_gene["start"] > e.start) & (ref_cds_gene["end"] == e.end)]
                    cds = cds.drop_duplicates(subset=["seq_id", "start", "end"], keep="first")
                    if cds.shape[0] == 1:
                        e.cds_s = cds.iloc[0, 1]
                        e.cds_e = cds.iloc[0, 2]
                elif strand == "-":
                    cds = ref_cds_gene[(ref_cds_gene["start"] == e.start) & (ref_cds_gene["end"] < e.end)]
                    cds = cds.drop_duplicates(subset=["seq_id", "start", "end"], keep="first")
                    if cds.shape[0] == 1:
                        e.cds_s = cds.iloc[0, 1]
                        e.cds_e = cds.iloc[0, 2]
                exon_dict[chr][gene][k][trans] = e

# annotation utr info using the reference
# +: seq_chr[k.start-1:k.start+2] == "ATG", seq_chr[k.end:k.end+3] in ["TGA", "TAA", "TAG"]
# -: seq_chr[k.end-3:k.end] == "CAT", seq_chr[k.start-4:k.start-1] in ["CTA", "TTA", "TCA"]
# "".join([str(seq_chr[x[0]-1:x[1]]) for x in exons])
ref_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(transcript_info)))
for index, row in ref_exon.iterrows():
    chr, start, end, strand, trans, gene, type = list(row)
    ref_dict[chr][gene+":"+strand][trans].type = type
    ref_dict[chr][gene+":"+strand][trans].strand = strand
    ref_dict[chr][gene+":"+strand][trans].exon.append(str(min(start, end))+"-"+str(max(start, end)))

for chr, i in ref_dict.items():
    for gene, j in i.items():
        strand = gene[-1]
        for trans, k in j.items():
            k.start, k.end = 0, 0
            for exon in k.exon:
                e = exon_dict[chr][gene][exon][trans]
                if e.cds_sm == 0:
                    continue
                if k.start == 0 or e.cds_sm < k.start:
                    k.start = e.cds_sm
                if k.end == 0 or e.cds_em > k.end:
                    k.end = e.cds_em
            if k.start == 0 and k.type != "protein_coding":
                k.start = int(k.exon[0].split("-")[0])
                k.end = int(k.exon[-1].split("-")[1])
            elif k.start == 0:
                print("???")
            ref_dict[chr][gene][trans] = k

trans_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for chr, i in ref_dict.items():
    for gene, j in i.items():
        strand = gene[-1]
        for trans, k in j.items():
            exons = [[int(y) for y in x.split("-")] for x in k.exon]
            len_iso = sum([x[1] - x[0] for x in exons])
            cds = []
            for e in exons:
                if e[1] <= k.start or e[0] >= k.end:
                    continue
                if e[0] <= k.start and e[1] >= k.end:
                    cds.append(f"{k.start}-{k.end}")
                    break
                if e[0] <= k.start and e[1] >= k.start:
                    cds.append(f"{k.start}-{e[1]}")
                elif e[0] <= k.end and e[1] >= k.end:
                    cds.append(f"{e[0]}-{k.end}")
                elif e[0] > k.start and e[1] < k.end:
                    cds.append(f"{e[0]}-{e[1]}")
            cds_info = ";".join(cds)
            trans_dict[chr][gene][cds_info].append(trans)
            len_utr = [0, 0]
            for es, ee in exons:
                if ee <= k.start:
                    len_utr[0] += ee - es
                else:
                    len_utr[0] += max(0, k.start - es)
                    break
            for es, ee in exons[::-1]:
                if k.end <= es:
                    len_utr[1] += ee - es
                else:
                    len_utr[1] += max(0, ee - k.end)
                    break
            if gene[-1] == "+":
                len_utr = len_utr[::-1]
            k.utr3 = len_utr[0]
            k.utr5 = len_utr[1]
            k.body = len_iso - len_utr[0] - len_utr[1]
            k.cds = cds_info
            k.id = ",".join(k.exon)
            if len(k.exon) > 1:
                k.id = "-".join(k.id.split("-")[1:-1])
            ref_dict[chr][gene][trans] = k

# def reads & novel
reads_dict = {}
undef_reads = {}
with open("/mnt/md0/oe_full_length/output/OEfulllength/isoquant_bam/OUT/OUT.read_assignments.tsv", "r") as fs:
    line = fs.readline()
    line = fs.readline()
    line = fs.readline()
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        if line[3] == "." or line[1] not in chr_list:
            line = fs.readline()
            continue
        if line[5] == "unique" and "mono_exon_match" in line[6]:
            t = ref_dict[line[1]][line[4]+":"+line[2]][line[3]]
            start, end = line[7].split(",")[0].split("-")
            if int(start) <= t.start and int(end) >= t.end:
                reads_dict[line[0]] = [line[1], line[2], line[4], line[3], line[7]]
                line = fs.readline()
                continue
        if line[4][:7] != "ENSMUSG" or len(line[7].split(",")) < 2:
            line = fs.readline()
            continue
        ts = line[7].split("-")
        ts_s = ts[0]
        ts_e = ts[-1]
        ts = "-".join(ts[1:-1])
        rec = ["", []]
        same_test = False
        gene = line[4]+":"+line[2]
        for k, v in ref_dict[line[1]][gene].items():
            if ts in v.id:
                if rec[0] != "":
                    if len(v.id) < len(rec[0]):
                        rec = [v.id, [k]]
                    elif v.id == rec[0]:
                        rec[1].append(k)
                    same_test = True
                else:
                    rec = [v.id, [k]]
        if same_test and ts != rec[0]:
            rec = ["", []]
        if rec[0] != "":
            #if len(rec[1]) > 1:
            #    print("???")
            #    break
            reads_dict[line[0]] = [line[1], line[2], line[4], trans_dict[line[1]][gene][ref_dict[line[1]][gene][rec[1][0]].cds][0], line[7]]
            line = fs.readline()
            continue
        if "tes_match" in line[6] and len(line[7].split(",")) > 1 and line[0] not in undef_reads:
            undef_reads[line[0]] = [line[1], line[4]+":"+line[2], line[7]]
        line = fs.readline()

undef_isoforms = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(str)))))
with open("/mnt/md0/oe_full_length/output/OEfulllength/read_tags.tsv", "r") as fs:
    line = fs.readline()
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        if line[0] in reads_dict or line[0] not in undef_reads:
            line = fs.readline()
            continue
        chr, gene, t = undef_reads[line[0]]
        t = t.split(",")
        t[0] = t[0].split("-")[1]
        t[-1] = t[-1].split("-")[0]
        undef_isoforms[chr][gene][",".join(t)][line[2]][line[5]] = line[0]
        line = fs.readline()

for chr, i in undef_isoforms.items():
    for gene, j in i.items():
        count = 0
        for term, k in j.items():
            if len(k) < 2 or max([len(x) for x in k.values()]) < 2:
                continue
            if (gene[-1] == "+" and term.split(",")[0] in [x.exon[0].split("-")[1] for x in ref_dict[chr][gene].values()]) or (gene[-1] == "-" and term.split(",")[-1] in [x.exon[-1].split("-")[0] for x in ref_dict[chr][gene].values()]):
                es = 0
                ee = 0
                for bc, x in k.items():
                    for umi, y in x.items():
                        e = [[int(y) for y in x.split("-")] for x in undef_reads[y][2].split(",")]
                        if es == 0 or e[0][0] > es:
                            es = e[0][0]
                        if ee == 0 or e[-1][1] < ee:
                            ee = e[-1][1]
                t = transcript_info()
                t.type = "NOVEL"
                t.exon = (str(es)+"-"+term+"-"+str(ee)).split(",")
                t.start = es
                t.end = ee
                exons = [[int(y) for y in x.split("-")] for x in t.exon]
                t.body = sum([x[1] - x[0] for x in exons])
                id = gene[:-2]+"_NOVEL"+str(count).zfill(2)
                ref_dict[chr][gene][id] = t
                count += 1
                if count > 100:
                    print("???")
                for bc, x in k.items():
                    for umi, y in x.items():
                        reads_dict[y] = [chr, gene[-1], gene[:-2], id, undef_reads[y][2]]

res = open("/mnt/md0/oe_full_length/output/OEfulllength/isoquant_bam/OUT/OUT.read_assignments_fix.tsv", "w")
with open("/mnt/md0/oe_full_length/output/OEfulllength/read_tags.tsv", "r") as fs:
    line = fs.readline()
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        if line[0] not in reads_dict:
            line = fs.readline()
            continue
        infos = reads_dict[line[0]]
        if infos[0] not in chr_list:
            line = fs.readline()
            continue
        trans_info = ref_dict[infos[0]][infos[2]+":"+infos[1]][infos[3]]
        exons = [[int(y) for y in x.split("-")] for x in infos[4].split(",")]
        len_iso = 0
        for es, ee in exons:
            len_iso += ee - es
        len_utr = [0, 0]
        for es, ee in exons:
            if ee <= trans_info.start:
                len_utr[0] += ee - es
            else:
                len_utr[0] += max(0, trans_info.start - es)
                break
        for es, ee in exons[::-1]:
            if trans_info.end <= es:
                len_utr[1] += ee - es
            else:
                len_utr[1] += max(0, ee - trans_info.end)
                break
        #if (infos[1] == "-" and len_utr[0] == 0) or (infos[1] == "+" and len_utr[1] == 0):
        #    print("???")
        #    break
        if infos[1] == "+":
            len_utr = len_utr[::-1]
        r = res.write(line[0]+"\t"+"\t".join(infos + [str(len_utr[0]), str(len_utr[1]), str(len_iso - len_utr[0] - len_utr[1]), line[2], line[5]])+"\n")
        line = fs.readline()

res.close()

cells = {}
cell_id = 0
with open("/mnt/md0/oe_full_length/output/OEfulllength/neuron_bc_info.tsv", "r") as fs:
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        cells[line[0]] = [cell_id, line[1]]
        cell_id += 1
        line = fs.readline()

trans_info = defaultdict(list)
trans = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
with open("/mnt/md0/oe_full_length/output/OEfulllength/isoquant_bam/OUT/OUT.read_assignments_fix.tsv", "r") as fs:
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        bc = line[9] + "-1"
        if bc not in cells:
            line = fs.readline()
            continue
        if line[4] not in trans_info:
            trans_info[line[4]] = [len(trans_info), line[1], line[2], line[3]]
        if line[10] not in trans[bc][line[4]]:
            trans[bc][line[4]][line[10]] = [int(line[6]), int(line[7]), int(line[8])]
        else:
            trans_len = trans[bc][line[4]][line[10]]
            trans[bc][line[4]][line[10]] = [max(int(line[6]), trans_len[0]), max(int(line[7]), trans_len[1]), max(int(line[8]), trans_len[2])]
        line = fs.readline()

cds_dict = {}
cds_info = defaultdict(list)
for chr, i in ref_dict.items():
    for gene, j in i.items():
        for t, k in j.items():
            if k.type != "protein_coding":
                continue
            cds = []
            for exon in k.exon:
                s, e = [int(x) for x in exon.split("-")]
                if k.start <= s and k.end >= e:
                    cds.append(f"{s}-{e}")
                elif k.start > s and k.end >= e:
                    cds.append(f"{k.start}-{e}")
                elif k.start <= s and k.end < e:
                    cds.append(f"{s}-{k.end}")
                elif k.start > s and k.end < e:
                    cds.append(f"{k.start}-{k.end}")
            group = "_".join([chr, gene, ";".join(cds)])
            cds_dict[t] = [chr, gene[-1], gene[:-2], group]
            cds_info[group].append(t)

count = [0, 0]
for chr, i in ref_dict.items():
    for gene, j in i.items():
        ass = {}
        ss = {}
        es = {}
        for t, k in j.items():
            if "_NOVEL" in t or len(k.exon) < 2:
                continue
            ss[k.exon[0].split("-")[1]] = 0
            es[k.exon[-1].split("-")[0]] = 0
            for x in ";".join(k.exon).split("-")[1:-1]:
                ass[x] = 0
        for t, k in j.items():
            if "_NOVEL" not in t or len(k.exon) < 2:
                continue
            if k.exon[0].split("-")[1] not in ss or k.exon[-1].split("-")[0] not in es:
                continue
            count[0] += 1
            ref_dict[chr][gene][t].type = "NOVEL_IC"
            for x in ";".join(k.exon).split("-")[1:-1]:
                if x not in ass:
                    ref_dict[chr][gene][t].type = "NOVEL_NC"
                    count[1] += 1
                    break


as_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(as_info)))
for chr, i in ref_dict.items():
    for gene, j in i.items():
        for t, k in j.items():
            if len(k.exon) < 2:
                continue
            terms = []
            for x in k.exon:
                terms += x.split("-")
            terms = [int(x) for x in set(terms[1:-1])]
            for x in terms:
                as_dict[chr][gene][x].trans.append(t)
                if k.utr3 > 0 and k.utr5 > 0:
                    if x >= k.start and x <= k.end:
                        if as_dict[chr][gene][x].type == "":
                            as_dict[chr][gene][x].type = "CDS"
                        elif as_dict[chr][gene][x].type == "UTR":
                            as_dict[chr][gene][x].type = "Hybrid"
                    else:
                        if as_dict[chr][gene][x].type == "":
                            as_dict[chr][gene][x].type = "UTR"
                        elif as_dict[chr][gene][x].type == "CDS":
                            as_dict[chr][gene][x].type = "Hybrid"

terms = {}
with open("/mnt/md0/oe_full_length/output/OEfulllength/read_tags.tsv", "r") as fs:
    line = fs.readline()
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        if line[0] in reads_dict:
            line = fs.readline()
            continue
        terms[line[0]] = [line[2], line[5]]
        line = fs.readline()

read_ass = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
with open("/mnt/md0/oe_full_length/output/OEfulllength/isoquant_bam/OUT/OUT.read_assignments.tsv", "r") as fs:
    line = fs.readline()
    line = fs.readline()
    line = fs.readline()
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        if line[3] == "." or line[1] not in chr_list:
            line = fs.readline()
            continue
        if line[0] not in terms or line[4][:7] != "ENSMUSG" or len(line[7].split(",")) < 2:
            line = fs.readline()
            continue
        gene = line[4]+":"+line[2]
        for x in [int(y) for y in set(line[7].replace(",", "-").split("-"))]:
            bc, umi = terms[line[0]]
            if x not in as_dict[line[1]][gene] or bc + "-1" not in cells:
                continue
            id = "_".join([line[1], gene, str(x)])
            read_ass[id][bc][umi] = 0
        line = fs.readline()

features = defaultdict(list)
feature_count = 0
count_row = []
count_col = []
count_data = []
for id, v in read_ass.items():
    chr, gene, strand, pos = id.replace(":", "_").split("_")
    features["id"].append(id)
    features["chr"].append(chr)
    features["strand"].append(strand)
    features["pos"].append(int(pos))
    features["gene"].append(gene)
    features["trans"].append(";".join(as_dict[chr][gene+":"+strand][int(pos)].trans))
    features["type"].append(as_dict[chr][gene+":"+strand][int(pos)].type)
    for bc, u in v.items():
        count_row.append(feature_count)
        count_col.append(cells[bc+"-1"][0])
        count_data.append(len(u))
    feature_count += 1

count_matrix = sparse.csc_matrix((np.array(count_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
with tables.open_file("/mnt/md0/oe_full_length/output/OEfulllength/oe_junction_info.h5", "w", title="") as fs:
    group = fs.create_group("/", "matrix", "Datasets of count")
    feature_group = fs.create_group(group, "features", "Genes and other features measured")
    for k, v in features.items():
        r = fs.create_array(feature_group, k, np.array(v)) 
    r = fs.create_array(group, "genes", np.array(features["id"]))
    r = fs.create_array(group, "barcodes", np.array(list(cells.keys())))
    r = fs.create_array(group, "data", count_matrix.data)
    r = fs.create_array(group, "indices", count_matrix.indices)
    r = fs.create_array(group, "indptr", count_matrix.indptr)
    r = fs.create_array(group, "shape", count_matrix.shape)


######################################################################################
cds = {}
features = defaultdict(list)
len_info = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for k, v in trans.items():
    for t, u in v.items():
        if t not in cds_dict:
            continue
        chr, strand, gene, group = cds_dict[t]
        if not group in cds:
            features["name"].append(group)
            features["id"].append(len(cds))
            features["chr"].append(chr)
            features["strand"].append(strand)
            features["gene"].append(gene)
            features["trans"].append(";".join(cds_info[group]))
            cds[group] = len(cds)
        if not k in len_info[k]:
            len_info[k][group]["count"] = [0]
        len_info[k][group]["count"][0] += len(u)
        for x in u.values():
            if x[0] > 0:
                len_info[k][group]["utr3"].append(x[0])
            if x[1] > 0:
                len_info[k][group]["utr5"].append(x[1])
            if x[2] > 0:
                len_info[k][group]["body"].append(x[2])

count_row = []
count_col = []
count_data = []
utr3_data = []
utr5_data = []
body_data = []
for k, v in len_info.items():
    for t, u in v.items():
        count_row.append(cds[t])
        count_col.append(cells[k][0])
        count_data.append(u["count"][0])
        if len(u["utr3"]) > 0:
            utr3_data.append(max(u["utr3"]))
        else:
            utr3_data.append(0)
        if len(u["utr5"]) > 0:
            utr5_data.append(max(u["utr5"]))
        else:
            utr5_data.append(0)
        if len(u["body"]) > 0:
            body_data.append(max(u["body"]))
        else:
            body_data.append(0)

count_matrix = sparse.csc_matrix((np.array(count_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
utr3_matrix = sparse.csc_matrix((np.array(utr3_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
utr5_matrix = sparse.csc_matrix((np.array(utr5_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
body_matrix = sparse.csc_matrix((np.array(body_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
with tables.open_file("/mnt/md0/oe_full_length/output/OEfulllength/osn_gene_info.h5", "w", title="") as fs:
    group = fs.create_group("/", "matrix", "Datasets of count")
    feature_group = fs.create_group(group, "features", "Genes and other features measured")
    for k, v in features.items():
        r = fs.create_array(feature_group, k, np.array(v))
    r = fs.create_array(group, "genes", np.array(features["name"]))
    r = fs.create_array(group, "barcodes", np.array(list(cells.keys())))
    r = fs.create_array(group, "data", count_matrix.data)
    r = fs.create_array(group, "utr3", utr3_matrix.data)
    r = fs.create_array(group, "utr5", utr5_matrix.data)
    r = fs.create_array(group, "body", body_matrix.data)
    r = fs.create_array(group, "indices", count_matrix.indices)
    r = fs.create_array(group, "indptr", count_matrix.indptr)
    r = fs.create_array(group, "shape", count_matrix.shape)


len_info = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for k, v in trans.items():
    for t, u in v.items():
        if t not in cds_dict:
            continue
        chr, strand, gene, group = cds_dict[t]
        utr3 = []
        utr5 = []
        body = []
        for x in u.values():
            if x[0] > 0:
                utr3.append(x[0])
            if x[1] > 0:
                utr5.append(x[1])
            if x[2] > 0:
                body.append(x[2])
        if len(utr3) > 0:
            len_info[cells[k][1]][group]["utr3"].append(max(utr3))
        if len(utr5) > 0:
            len_info[cells[k][1]][group]["utr5"].append(max(utr5))
        if len(body) > 0:
            len_info[cells[k][1]][group]["body"].append(max(body))
        if len(utr3) > 0 and len(utr5) > 0 and len(body) > 0:
            len_info[cells[k][1]][group]["len"].append(max(utr3)+max(utr5)+max(body))

with open("/mnt/md0/oe_full_length/output/OEfulllength/osn_gene_len.tsv", "w") as fs:
    r = fs.write("\t".join(["type", "trans", "utr3", "utr5", "body", "len", "r3", "r5", "rlen"])+"\n")
    for k, v in len_info.items():
        for t, u in v.items():
            utr3 = round(sum(u["utr3"])/len(u["utr3"])) if len(u["utr3"]) > 0 else 0
            utr5 = round(sum(u["utr5"])/len(u["utr5"])) if len(u["utr5"]) > 0 else 0
            body = round(sum(u["body"])/len(u["body"])) if len(u["body"]) > 0 else 0
            length = round(sum(u["len"])/len(u["len"])) if len(u["len"]) > 0 else 0
            #tif = trans_info[t]
            #tif = ref_dict[tif[1]][tif[3]+":"+tif[2]][t]
            r = fs.write("\t".join([k, t, str(utr3), str(utr5), str(body), str(length), "0", "0", "0"])+"\n")

features = defaultdict(list)
for k, v in trans_info.items():
    features["name"].append(k)
    features["id"].append(v[0])
    features["chr"].append(v[1])
    features["strand"].append(v[2])
    features["gene"].append(v[3])
    features["type"].append(ref_dict[v[1]][v[3]+":"+v[2]][k].type)
    features["utr3"].append(ref_dict[v[1]][v[3]+":"+v[2]][k].utr3)
    features["utr5"].append(ref_dict[v[1]][v[3]+":"+v[2]][k].utr5)
    features["body"].append(ref_dict[v[1]][v[3]+":"+v[2]][k].body)
    start = ref_dict[v[1]][v[3]+":"+v[2]][k].start
    end = ref_dict[v[1]][v[3]+":"+v[2]][k].end
    exons = [[int(y) for y in x.split("-")] for x in ref_dict[v[1]][v[3]+":"+v[2]][k].exon]
    exon_rec = [[], [], []]
    for s, e in exons:
        if e <= start:
            exon_rec[0].append(f"{s}-{e}")
        elif s >= end:
            exon_rec[2].append(f"{s}-{e}")
        elif s < start and end < e:
            exon_rec[0].append(f"{s}-{start}")
            exon_rec[1].append(f"{start}-{end}")
            exon_rec[2].append(f"{end}-{e}")
        elif start > s and start < e:
            exon_rec[0].append(f"{s}-{start}")
            exon_rec[1].append(f"{start}-{e}")
        elif s < end and e > end:
            exon_rec[1].append(f"{s}-{end}")
            exon_rec[2].append(f"{end}-{e}")
        else:
            exon_rec[1].append(f"{s}-{e}")
    features["exon"].append(";".join([",".join(x) for x in exon_rec]))

count_row = []
count_col = []
count_data = []
utr3_data = []
utr5_data = []
for k, v in trans.items():
    for t, u in v.items():
        count_row.append(trans_info[t][0])
        count_col.append(cells[k][0])
        count_data.append(len(u))
        utr3 = []
        utr5 = []
        for umi, x in u.items():
            if x[0] > 0:
                utr3.append(x[0])
            if x[1] > 0:
                utr5.append(x[1])
        if len(utr3) > 0:
            utr3_data.append(round(np.median(utr3)))
        else:
            utr3_data.append(0)
        if len(utr5) > 0:
            utr5_data.append(round(np.median(utr5)))
        else:
            utr5_data.append(0)

count_matrix = sparse.csc_matrix((np.array(count_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
utr3_matrix = sparse.csc_matrix((np.array(utr3_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
utr5_matrix = sparse.csc_matrix((np.array(utr5_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
with tables.open_file("/mnt/md0/oe_full_length/output/OEfulllength/osn_trans_ass.h5", "w", title="") as fs:
    group = fs.create_group("/", "matrix", "Datasets of count")
    feature_group = fs.create_group(group, "features", "Genes and other features measured")
    for k, v in features.items():
        r = fs.create_array(feature_group, k, np.array(v))
    r = fs.create_array(group, "genes", np.array(features["name"]))
    r = fs.create_array(group, "barcodes", np.array(list(cells.keys())))
    r = fs.create_array(group, "data", count_matrix.data)
    r = fs.create_array(group, "utr3", utr3_matrix.data)
    r = fs.create_array(group, "utr5", utr5_matrix.data)
    r = fs.create_array(group, "indices", count_matrix.indices)
    r = fs.create_array(group, "indptr", count_matrix.indptr)
    r = fs.create_array(group, "shape", count_matrix.shape)

len_info = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for k, v in trans.items():
    for t, u in v.items():
        for x in u.values():
            if x[0] > 0:
                len_info[cells[k][1]][t]["utr3"].append(x[0])
            if x[1] > 0:
                len_info[cells[k][1]][t]["utr5"].append(x[1])
            if x[2] > 0:
                len_info[cells[k][1]][t]["body"].append(x[2])

with open("/mnt/md0/oe_full_length/output/OEfulllength/osn_iso_len.tsv", "w") as fs:
    r = fs.write("\t".join(["type", "trans", "utr3", "utr5", "len", "r3", "r5", "rlen"])+"\n")
    for k, v in len_info.items():
        for t, u in v.items():
            utr3 = round(sum(u["utr3"])/len(u["utr3"])) if len(u["utr3"]) > 0 else 0
            utr5 = round(sum(u["utr5"])/len(u["utr5"])) if len(u["utr5"]) > 0 else 0
            body = round(sum(u["body"])/len(u["body"])) if len(u["body"]) > 0 else 0
            tif = trans_info[t]
            tif = ref_dict[tif[1]][tif[3]+":"+tif[2]][t]
            r = fs.write("\t".join([k, t, str(utr3), str(utr5), str(body+utr3+utr5), str(tif.utr3), str(tif.utr5), str(tif.body + tif.utr3 + tif.utr5)])+"\n")

with open("/mnt/md0/oe_full_length/output/OEfulllength/osn_trans_len.tsv", "w") as fs:
    r = fs.write("\t".join(["type", "trans", "utr3", "utr5", "len", "r3", "r5", "rlen"])+"\n")
    for k, v in trans.items():
        for t, u in v.items():
            for x in u.values():
                tif = trans_info[t]
                tif = ref_dict[tif[1]][tif[3]+":"+tif[2]][t]
                r = fs.write("\t".join([cells[k][1], t, str(x[0]), str(x[1]), str(sum(x)), str(tif.utr3), str(tif.utr5), str(tif.body + tif.utr3 + tif.utr5)])+"\n")

gene_iso = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
for k, v in trans.items():
    for t, u in v.items():
        for x in u.values():
            tif = trans_info[t]
            if ref_dict[tif[1]][tif[3]+":"+tif[2]][t].type == "protein_coding":
                gene_iso[cells[k][1]][tif[3]+":"+tif[2]][t] = 0

with open("/mnt/md0/oe_full_length/output/OEfulllength/osn_gene_iso.tsv", "w") as fs:
    r = fs.write("\t".join(["type", "gene", "count"])+"\n")
    for k, v in gene_iso.items():
        for g, t in v.items():
            c = 0
            for s in t.keys():
                info = trans_info[s]
                if ref_dict[info[1]][g][s].type == "protein_coding":
                    c += 1
            r = fs.write("\t".join([k, g, str(c)])+"\n")


samfile = "/mnt/md0/oe_full_length/output/OEfulllength/tagged.bam"
chr_list = ["chr"+str(x) for x in list(range(1, 50))+["X", "Y"]]
ref = pd.read_table("~/library/refdata-gex-mm10-2020-A/genes/genes.gtf", comment="#", header=None, dtype={0:str})
ref.columns = ["seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
ref = ref[(ref["type"] == "gene") & (ref["seq_id"].isin(chr_list))]
ref["symbol"] = list(map(lambda x: re.findall(r'gene_name "(.*?)";', x)[0], ref["attributes"]))
ref["attributes"] = list(map(lambda x: re.findall(r'gene_id "(.*?)";', x)[0], ref["attributes"]))
ref = ref[["seq_id", "start", "end", "strand", "symbol", "attributes"]].sort_values(by=["seq_id", "start", "end", "attributes"])
ref.drop_duplicates(subset=["seq_id", "start", "end", "attributes"], keep="first", inplace=True)

cells = {}
cell_id = 0
cell_types = defaultdict(int)
with open("/mnt/md0/oe_full_length/output/OEfulllength/neuron_bc_info.tsv", "r") as fs:
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        cells[line[0]] = [cell_id, line[1]]
        cell_id += 1
        cell_types[line[1]] += 1
        line = fs.readline()

terms = {}
with open("/mnt/md0/oe_full_length/output/OEfulllength/isoquant_bam/OUT/OUT.read_assignments_fix.tsv", "r") as fs:
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        bc = line[9] + "-1"
        if bc not in cells:
            line = fs.readline()
            continue
        term = "_".join([line[4], line[9], line[10]])
        term_len = int(line[6]) + int(line[7]) + int(line[8])
        if term in terms:
            if term_len > terms[term][5]:
                terms[term] = line[:5]+[term_len]
        else:
            terms[term] = line[:5]+[term_len]
        line = fs.readline()

term_info = {}
for k, v in terms.items():
    term_info[v[0]] = v[1:]

with pysam.AlignmentFile(samfile) as fs:
    res = pysam.AlignmentFile("/mnt/md0/oe_full_length/output/OEfulllength/tagged_fix.bam", "wb", header=fs.header)
    for reads in fs.fetch():
        if reads.qname not in term_info:
            continue
        info = term_info[reads.qname]
        reads.set_tag("GN", info[2])
        reads.set_tag("TR", info[3])
        r = res.write(reads)
    res.close()

r = os.system("samtools sort -@8 tagged_fix.bam -o tagged_fix.sort.bam")
r = os.system("samtools index tagged_fix.sort.bam")
r = os.system("rm tagged_fix.bam")
r = os.system("scp tagged_fix.sort.bam* xuyc7@202.116.90.56:/usr/share/nginx/html/jbrowse2/usr_data/")
with pysam.AlignmentFile("/mnt/md0/oe_full_length/output/OEfulllength/tagged_fix.sort.bam") as fs:
    res_dict = {}
    for type in cell_types.keys():
        res_dict[type] = pysam.AlignmentFile(f"/mnt/md0/oe_full_length/output/OEfulllength/tagged_fix_{type}.bam", "wb", header=fs.header)
    for reads in fs.fetch():
        bc = reads.get_tag("CB") + "-1"
        if bc in cells:
            r = res_dict[cells[bc][1]].write(reads)
    for type in cell_types.keys():
        res_dict[type].close()
        r = os.system(f"samtools index /mnt/md0/oe_full_length/output/OEfulllength/tagged_fix_{type}.bam")


ref_sub = ref[ref["symbol"].str.contains(r"^Olfr", regex=True)]
#/mnt/md0/oe_full_length/output/OEfulllength/olfr_info.csv
olfr_res = open("olfr_info.csv", "w")
r = olfr_res.write(f"ID,Gene,Symbol,TypeC,TypeE,TSS,UTR,CDS,Exon\n")
for index, row in ref_sub.iterrows():
    chr, start, end, strand, symbol, gene_id = list(row)
    gene_id += ":" + strand
    if len(ref_dict[chr][gene_id]) < 2:
        t, k = list(ref_dict[chr][gene_id].items())[0]
        cds = [0, 0]
        for s, e in [[int(y) for y in x.split("-")] for x in k.cds.split(";")]:
            if e - s > cds[1] - cds[0]:
                cds = [s, e]
        for s, e in [[int(y) for y in x.split("-")] for x in k.cds.split(";")]:
            if [s, e] != cds and e - s > (cds[1] - cds[0]) * 0.4:
                print("???")
        r = olfr_res.write(f"{t},{gene_id},{symbol},-,-,-,-,{k.cds},{';'.join(k.exon)}\n")
        continue
    term_main = ["", 0, 0, 0, 0, 0]
    for t, k in ref_dict[chr][gene_id].items():
        if k.utr5 == 0 or k.utr3 == 0 or (term_main[0] != "" and k.body < term_main[1].body):
            continue
        term_main = [t, k, 0, 0, 0, 0]
    if term_main[0] == "":
        for t, k in ref_dict[chr][gene_id].items():
            r = olfr_res.write(f"{t},{gene_id},{symbol},Unknown,Unknown,Unknown,Unknown,{k.cds},{';'.join(k.exon)}\n")
        continue
    cds = [0, 0]
    for s, e in [[int(y) for y in x.split("-")] for x in term_main[1].cds.split(";")]:
        if e - s > cds[1] - cds[0]:
            cds = [s, e]
    for t, k in ref_dict[chr][gene_id].items():
        if k.utr5 == 0 or k.utr3 == 0:
            continue
        for es, ee in [[int(y) for y in x.split("-")] for x in k.exon]:
            if es <= cds[0] and ee >= cds[1] and ee - es > term_main[5] - term_main[4]:
                term_main = [t, k, cds[0], cds[1], es, ee]
    tsss, tsse = [int(x) for x in term_main[1].exon[0].split("-")] if strand == "+" else [int(x) for x in term_main[1].exon[-1].split("-")]
    utrs, utre = [], []
    for exon in term_main[1].exon:
        s, e = [int(x) for x in exon.split("-")]
        if s < term_main[1].start:
            utrs.append(f"{s}-{min(term_main[1].start, e)}")
        if e > term_main[1].end:
            utre.append(f"{max(term_main[1].end, s)}-{e}")
    #mutr3 = "-".join(";".join(utrs).split("-")[1:-1]) if strand == "-" else "-".join(";".join(utre).split("-")[1:-1])
    #mutr5 = "-".join(";".join(utrs).split("-")[1:-1]) if strand == "+" else "-".join(";".join(utre).split("-")[1:-1])
    mutr3 = ";".join(utrs) if strand == "-" else ";".join(utre)
    mutr5 = ";".join(utrs) if strand == "+" else ";".join(utre)
    for t, k in ref_dict[chr][gene_id].items():
        if t == term_main[0]:
            r = olfr_res.write(f"{t},{gene_id},{symbol},-,-,-,-,{k.cds},{';'.join(k.exon)}\n")
            continue
        type = ["", "", "-", "-"]
        # TSS
        if strand == "+":
            s, e = [int(x) for x in k.exon[0].split("-")]
            if e <= tsss or s > tsse:
                type[2] = "ATSS"
        else:
            s, e = [int(x) for x in k.exon[-1].split("-")]
            if e <= tsss or s > tsse:
                type[2] = "ATSS"
        # UTR
        utrs, utre = [], []
        for exon in k.exon:
            s, e = [int(x) for x in exon.split("-")]
            if s < term_main[1].start:
                utrs.append(f"{s}-{min(term_main[1].start, e)}")
            if e > term_main[1].end:
                utre.append(f"{max(term_main[1].end, s)}-{e}")
        #utr3 = "-".join(";".join(utrs).split("-")[1:-1]) if strand == "-" else "-".join(";".join(utre).split("-")[1:-1])
        #utr5 = "-".join(";".join(utrs).split("-")[1:-1]) if strand == "+" else "-".join(";".join(utre).split("-")[1:-1])
        utr3 = ";".join(utrs) if strand == "-" else ";".join(utre)
        utr5 = ";".join(utrs) if strand == "+" else ";".join(utre)
        if utr3 != mutr3 and utr5 != mutr5:
            type[3] = "UTR3+UTR5"
        elif utr3 != mutr3:
            type[3] = "UTR3"
        elif utr5 != mutr5:
            type[3] = "UTR5"
        # CDS: Full
        for exon in k.exon:
            s, e = [int(x) for x in exon.split("-")]
            if s <= term_main[2] and term_main[3] <= e:
                type[0] = "-"
                break
        # CDS: SE
        if type[0] == "":
            for exon in k.exon:
                s, e = [int(x) for x in exon.split("-")]
                if (term_main[2] <= s and term_main[3] >= s) or (term_main[2] <= e and term_main[3] >= e):
                    type[0] = "X"
                    break
            type[0] = "" if type[0] == "X" else "SE"
        # CDS: RI, A3SS, A5SS
        if type[0] == "":
            cs = 0
            for exon in k.exon:
                s, e = [int(x) for x in exon.split("-")]
                if s <= term_main[2] and e >= term_main[2]:
                    type[0] += "A"
                if s <= term_main[3] and e >= term_main[3]:
                    type[0] += "B"
                if term_main[2] <= s and term_main[3] >= e:
                    cs += 1
            if type[0] == "AB" or type[0] == "BA":
                type[0] = "RI"
            elif (type[0] == "A" and strand == "+") or (type[0] == "B" and strand == "-"):
                type[0] = "A5SS" if cs == 0 else "A5SS+RI"
            elif (type[0] == "B" and strand == "+") or (type[0] == "A" and strand == "-"):
                type[0] = "A3SS" if cs == 0 else "A3SS+RI"
            elif type[0] == "":
                type[0] = "A3SS+A5SS" if cs == 1 else "A3SS+RI"
        # Exon: Full
        for exon in k.exon:
            s, e = [int(x) for x in exon.split("-")]
            if s <= term_main[4] and term_main[5] <= e:
                type[1] = "-"
                break
        # Exon: SE
        if type[1] == "":
            for exon in k.exon:
                s, e = [int(x) for x in exon.split("-")]
                if (term_main[4] <= s and term_main[5] >= s) or (term_main[4] <= e and term_main[5] >= e):
                    type[1] = "X"
                    break
            type[1] = "" if type[1] == "X" else "SE"
        # Exon: RI, A3SS, A5SS
        if type[1] == "":
            cs = 0
            for exon in k.exon:
                s, e = [int(x) for x in exon.split("-")]
                if s <= term_main[4] and e >= term_main[4]:
                    type[1] += "A"
                if s <= term_main[5] and e >= term_main[5]:
                    type[1] += "B"
                if term_main[4] <= s and term_main[5] >= e:
                    cs += 1
            if type[1] == "AB" or type[1] == "BA":
                type[1] = "RI"
            elif (type[1] == "A" and strand == "+") or (type[1] == "B" and strand == "-"):
                type[1] = "A5SS" if cs == 0 else "A5SS+RI"
            elif (type[1] == "B" and strand == "+") or (type[1] == "A" and strand == "-"):
                type[1] = "A3SS" if cs == 0 else "A3SS+RI"
            elif type[1] == "":
                type[1] = "A3SS+A5SS" if cs == 1 else "A3SS+RI"
        if type[0] == "" or type[1] == "":
            print(t, chr, symbol, gene_id)
        r = olfr_res.write(f"{t},{gene_id},{symbol},{type[0]},{type[1]},{type[2]},{type[3]},{k.cds},{';'.join(k.exon)}\n")

olfr_res.close()


def dna2protein(seq):
    codon_table = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M','ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T','AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K','AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R','CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L','CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P','CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q','CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R','GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V','GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A','GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E','GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G','TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S','TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L','TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_','TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
    protein = ""
    if len(seq) % 3 != 0:
        return ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3].upper()
        amino_acid = codon_table.get(codon, '?')
        if amino_acid == '_':
            break
        protein += amino_acid
    return protein

res = open("/mnt/md0/oe_full_length/output/OEfulllength/iso_trans_test.fa", "w")
res_aa = open("/mnt/md0/oe_full_length/output/OEfulllength/iso_trans_test.faa", "w")
iso_trans = pd.read_table("/home/cbh/work/oe_fl_re/oe_iso_trans_test.tsv")
for index, row in iso_trans.iterrows():
    chr = row["Chr"]
    gene = row["Gene"]
    symbol = row["Symbol"]
    id = row["ID"]
    strand = row["Strand"]
    seq = ""
    for s, e in [[int(y) for y in x.split("-")] for x in row["CDS"].split(";")]:
        seq += str(genome[chr].seq[s-1:e])
    if strand == "-":
        seq = seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdN', 'TGCAtgcaYRKMyrkmBVDHbvdhN'))[::-1]
    seq = seq[:len(seq)//3*3]
    seq_aa = dna2protein(seq)
    r = res.write(f">{symbol}_{chr}_{gene}_{id}\n{seq}\n")
    r = res_aa.write(f">{symbol}_{chr}_{gene}_{id}\n{seq_aa}\n")

res.close()
res_aa.close()
r = os.system("clustalo -i iso_trans_test.fa -o iso_trans_test.aln --outfmt=clu --distmat-out=iso_trans_test.mat --full --force")
r = os.system("clustalo -i iso_trans_test.faa -o iso_trans_test_aa.aln --outfmt=clu --distmat-out=iso_trans_test_aa.mat --full --force")


terms = defaultdict(str)
with open("iso_trans_test_aa.aln", "r") as fs:
    line = fs.readline()
    line = fs.readline()
    line = fs.readline()
    line = fs.readline()
    while line:
        line_raw = line
        line = [x.strip() for x in line.strip().split(" ")]
        if len(line) < 2 or min(len(line[0]), len(line[-1])) == 0 or line[0] == "*":
            line = fs.readline()
            continue
        terms[line[0]] += line[-1]
        line = fs.readline()

with open("iso_trans_test_aa_cal.tsv", "w") as fs:
    r = fs.write("\t".join(["Term", "Domain1", "Domain2", "Domain3", "Domain4", "Domain5", "Domain6", "Domain7"])+"\n")
    for k, v in terms.items():
        t = [k]
        t.append(str(round(100-v[29:58].count("-")*100/29, 2)))
        t.append(str(round(100-v[64:94].count("-")*100/30, 2)))
        t.append(str(round(100-v[104:138].count("-")*100/34, 2)))
        t.append(str(round(100-v[151:176].count("-")*100/25, 2)))
        t.append(str(round(100-v[198:234].count("-")*100/36, 2)))
        t.append(str(round(100-v[257:278].count("-")*100/21, 2)))
        t.append(str(round(100-v[288:311].count("-")*100/23, 2)))
        r = fs.write("\t".join(t)+"\n")


cells = {}
with open("/home/cbh/work/oe_fl_re/oe_olfr_info.tsv", "r") as fs:
    line = fs.readline()
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        cells[line[0][:-2]] = line[2]
        line = fs.readline()

res_dict = {}
with open("/mnt/md0/oe_full_length/output/OEfulllength/isoquant_bam/OUT/OUT.read_assignments_fix.tsv", "r") as fs:
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        res_dict[line[0]] = 0
        line = fs.readline()

olfr_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
fs = pysam.AlignmentFile("/mnt/md0/oe_full_length/output/OEfulllength/tagged_fix.sort.bam")
res_fn = "/mnt/md0/oe_full_length/output/OEfulllength/tagged_olfr.bam"
res = pysam.AlignmentFile(res_fn, "wb", header=fs.header)
rec_record = {}
ref_sub = ref[ref["symbol"].str.contains(r"^Olfr", regex=True)]
for index, row in ref_sub.iterrows():
    chr, start, end, strand, symbol, gene_id = list(row)
    gene_id += ":" + strand
    for reads in fs.fetch(chr, start, end):
        if not reads.qname in res_dict:
            continue
        if reads.qname in rec_record:
            continue
        if not reads.get_tag("CB") in cells:
            continue
        r = res.write(reads)
        olfr_dict[gene_id][reads.get_tag("CB")][reads.get_tag("UB")] = 0
        rec_record[reads.qname] = 0

res.close()
fs.close()
r = os.system("samtools sort -@8 tagged_olfr.bam -o tagged_olfr.sort.bam")
r = os.system("samtools index tagged_olfr.sort.bam")
r = os.system("rm tagged_olfr.bam")
r = os.system("scp tagged_olfr.sort.bam* xuyc7@202.116.90.56:/usr/share/nginx/html/jbrowse2/usr_data/")
olfr_res = open("olfr_count.csv", "w")
r = olfr_res.write("Gene,ID,Pos,Count,PerCell\n")
for gene, terms in olfr_dict.items():
    count = 0
    for k, v in terms.items():
        count += len(v)
    chr, start, end, strand, symbol, id = list(ref[(ref["attributes"] == gene[:-2])].iloc[0])
    r = olfr_res.write(f"{symbol},{gene[:-2]},{chr}:{start}-{end},{count},{round(count/len(terms))}\n")

olfr_res.close()




###############################################################################################

trans_info = defaultdict(list)
trans = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
with open("/mnt/md0/oe_full_length/output/OEfulllength/isoquant_bam/OUT/OUT.read_assignments_fix.tsv", "r") as fs:
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        bc = line[9] + "-1"
        if bc not in cells:
            line = fs.readline()
            continue
        if line[3] not in trans_info:
            trans_info[line[3]] = [len(trans_info), line[1], line[2], line[3]]
        if line[10] not in trans[bc][line[3]]:
            trans[bc][line[3]][line[10]] = [int(line[6]), int(line[7]), int(line[8])]
        else:
            trans_len = trans[bc][line[3]][line[10]]
            trans[bc][line[3]][line[10]] = [max(int(line[6]), trans_len[0]), max(int(line[7]), trans_len[1]), max(int(line[8]), trans_len[2])]
        line = fs.readline()

features = defaultdict(list)
for k, v in trans_info.items():
    features["name"].append(k)
    features["id"].append(v[0])
    features["chr"].append(v[1])
    features["strand"].append(v[2])
    features["gene"].append(v[3])

count_row = []
count_col = []
count_data = []
utr3_data = []
utr5_data = []
for k, v in trans.items():
    for t, u in v.items():
        count_row.append(trans_info[t][0])
        count_col.append(cells[k][0])
        count_data.append(len(u))
        utr3 = []
        utr5 = []
        for umi, x in u.items():
            if x[0] > 0:
                utr3.append(x[0])
            if x[1] > 0:
                utr5.append(x[1])
        if len(utr3) > 0:
            utr3_data.append(round(np.median(utr3)))
        else:
            utr3_data.append(0)
        if len(utr5) > 0:
            utr5_data.append(round(np.median(utr5)))
        else:
            utr5_data.append(0)

count_matrix = sparse.csc_matrix((np.array(count_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
utr3_matrix = sparse.csc_matrix((np.array(utr3_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
utr5_matrix = sparse.csc_matrix((np.array(utr5_data), (np.array(count_row), np.array(count_col))), shape=(len(features["id"]), len(cells)))
with tables.open_file("/mnt/md0/oe_full_length/output/OEfulllength/oe_gene_ass.h5", "w", title="") as fs:
    group = fs.create_group("/", "matrix", "Datasets of count")
    feature_group = fs.create_group(group, "features", "Genes and other features measured")
    for k, v in features.items():
        r = fs.create_array(feature_group, k, np.array(v))
    r = fs.create_array(group, "genes", np.array(features["name"]))
    r = fs.create_array(group, "barcodes", np.array(list(cells.keys())))
    r = fs.create_array(group, "data", count_matrix.data)
    r = fs.create_array(group, "utr3", utr3_matrix.data)
    r = fs.create_array(group, "utr5", utr5_matrix.data)
    r = fs.create_array(group, "indices", count_matrix.indices)
    r = fs.create_array(group, "indptr", count_matrix.indptr)
    r = fs.create_array(group, "shape", count_matrix.shape)


types = defaultdict(int)
cells = {}
cell_id = 0
with open("/home/cbh/work/oe_full_length/neuron_bc_info.tsv", "r") as fs:
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        cells[line[0]] = [cell_id, line[1]]
        types[line[1]] += 1
        cell_id += 1
        line = fs.readline()

types = list(types.keys())
res_dict = {}
with open("/mnt/md0/oe_full_length/output/OEfulllength/isoquant_bam/OUT/OUT.read_assignments_fix.tsv", "r") as fs:
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        if line[9]+"-1" in cells:
            res_dict[line[0]] = cells[line[9]+"-1"][1]
        line = fs.readline()

fs = pysam.AlignmentFile("/mnt/md0/oe_full_length/output/OEfulllength/tagged.bam")
bam_list = {}
for x in types:
    bam_list[x] = pysam.AlignmentFile(f"/mnt/md0/oe_full_length/output/OEfulllength/tagged_{x}.bam", "wb", header=fs.header)

for reads in fs.fetch(chr, start, end):


for reads in fs.fetch():
    if reads.qname in res_dict:
        r = bam_list[res_dict[reads.qname]].write(reads)

fs.close()
for x in types:
    bam_list[x].close()
    r = os.system(f"samtools sort -@8 tagged_{x}.bam -o tagged_{x}.sort.bam")
    r = os.system(f"samtools index tagged_{x}.sort.bam")
    r = os.system(f"rm tagged_{x}.bam")

for x in types:
    r = os.system(f"scp tagged_{x}.sort.bam xuyc7@202.116.90.56:/usr/share/nginx/html/jbrowse2/usr_data/")
    r = os.system(f"scp tagged_{x}.sort.bam.bai xuyc7@202.116.90.56:/usr/share/nginx/html/jbrowse2/usr_data/")



chr_list = ["chr"+str(x) for x in list(range(1, 50))+["X", "Y"]]
ref = pd.read_table("~/library/refdata-gex-mm10-2020-A/genes/genes.gtf", comment="#", header=None, dtype={0:str})
ref.columns = ["seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
ref = ref[(ref["type"] == "gene") & (ref["seq_id"].isin(chr_list))]
ref["symbol"] = list(map(lambda x: re.findall(r'gene_name "(.*?)";', x)[0], ref["attributes"]))
ref["attributes"] = list(map(lambda x: re.findall(r'gene_id "(.*?)";', x)[0], ref["attributes"]))
ref = ref[["seq_id", "start", "end", "strand", "symbol", "attributes"]].sort_values(by=["seq_id", "start", "end", "attributes"])
ref.drop_duplicates(subset=["seq_id", "start", "end", "attributes"], keep="first", inplace=True)

types = defaultdict(int)
cells = {}
cell_id = 0
with open("/home/cbh/work/oe_full_length/neuron_bc_info.tsv", "r") as fs:
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        cells[line[0][:-2]] = [cell_id, line[1]]
        types[line[1]] += 1
        cell_id += 1
        line = fs.readline()

fs = pysam.AlignmentFile("/mnt/md0/oe_full_length/output/OEfulllength/tagged.bam")
res = pysam.AlignmentFile(f"/mnt/md0/oe_full_length/output/OEfulllength/tagged_neuron_hdac9_raw.bam", "wb", header=fs.header)
chr, start, end = list(ref[ref["symbol"] == "Hdac9"].iloc[0])[:3]
count = 0
for reads in fs.fetch(chr, start, end):
    if not reads.has_tag("CB"):
        continue
    cb = reads.get_tag("CB")
    if not cb in cells:
        continue
    if reads.qlen < 500:
        continue
    lr = 0
    for cigar in reads.cigar:
        if cigar[0] == 0 or cigar[0] == 3:
            lr += cigar[1]
    if lr < 10000:
        continue
    r = res.write(reads)
    count += 1

fs.close()
res.close()
r = os.system(f"samtools sort -@8 tagged_neuron_hdac9_raw.bam -o tagged_neuron_hdac9_raw.sort.bam")
r = os.system(f"samtools index tagged_neuron_hdac9_raw.sort.bam")
r = os.system(f"rm tagged_neuron_hdac9_raw.bam")

r = os.system(f"scp tagged_neuron_hdac9_raw.sort.bam xuyc7@202.116.90.56:/usr/share/nginx/html/jbrowse2/usr_data/")
r = os.system(f"scp tagged_neuron_hdac9_raw.sort.bam.bai xuyc7@202.116.90.56:/usr/share/nginx/html/jbrowse2/usr_data/")





cells_d = {}
cells_v = {}
with open("/home/cbh/work/oe_full_length/barcode_osn_dv.csv", "r") as fs:
    line = fs.readline()
    line = fs.readline()
    while line:
        line = line.strip().split(",")
        if line[1] == "D":
            cells_d[line[0][:-2]] = 0
        else:
            cells_v[line[0][:-2]] = 0
        line = fs.readline()

###############################################################################################


fs = pysam.AlignmentFile("/mnt/md0/oe_full_length/output/OEfulllength/tagged.bam")
res_d = pysam.AlignmentFile(f"/mnt/md0/oe_full_length/output/OEfulllength/tagged_neuron_hdac9_d_raw.bam", "wb", header=fs.header)
res_v = pysam.AlignmentFile(f"/mnt/md0/oe_full_length/output/OEfulllength/tagged_neuron_hdac9_v_raw.bam", "wb", header=fs.header)
chr, start, end = list(ref[ref["symbol"] == "Hdac9"].iloc[0])[:3]
count = 0
for reads in fs.fetch(chr, start, end):
    if not reads.has_tag("CB"):
        continue
    cb = reads.get_tag("CB")
    if cb not in cells_d and cb not in cells_v:
        continue
    if reads.qlen < 500:
        continue
    lr = 0
    for cigar in reads.cigar:
        if cigar[0] == 0 or cigar[0] == 3:
            lr += cigar[1]
    if lr < 10000:
        continue
    if cb in cells_d:
        r = res_d.write(reads)
    else:
        r = res_v.write(reads)

fs.close()
res_d.close()
res_v.close()
r = os.system(f"samtools sort -@8 tagged_neuron_hdac9_d_raw.bam -o tagged_neuron_hdac9_d_raw.sort.bam")
r = os.system(f"samtools index tagged_neuron_hdac9_d_raw.sort.bam")
r = os.system(f"rm tagged_neuron_hdac9_d_raw.bam")
r = os.system(f"scp tagged_neuron_hdac9_d_raw.sort.bam xuyc7@202.116.90.56:/usr/share/nginx/html/jbrowse2/usr_data/")
r = os.system(f"scp tagged_neuron_hdac9_d_raw.sort.bam.bai xuyc7@202.116.90.56:/usr/share/nginx/html/jbrowse2/usr_data/")

r = os.system(f"samtools sort -@8 tagged_neuron_hdac9_v_raw.bam -o tagged_neuron_hdac9_v_raw.sort.bam")
r = os.system(f"samtools index tagged_neuron_hdac9_v_raw.sort.bam")
r = os.system(f"rm tagged_neuron_hdac9_v_raw.bam")
r = os.system(f"scp tagged_neuron_hdac9_v_raw.sort.bam xuyc7@202.116.90.56:/usr/share/nginx/html/jbrowse2/usr_data/")
r = os.system(f"scp tagged_neuron_hdac9_v_raw.sort.bam.bai xuyc7@202.116.90.56:/usr/share/nginx/html/jbrowse2/usr_data/")





#samtools sort -@8 tagged_fix.bam -o tagged_fix.sort.bam
#samtools index tagged_fix.sort.bam
#scp tagged_fix.sort.bam xuyc7@202.116.90.56:/usr/share/nginx/html/jbrowse2/usr_data/
#scp tagged_fix.sort.bam.bai xuyc7@202.116.90.56:/usr/share/nginx/html/jbrowse2/usr_data/









res = [0, 0, 0]
with open("/mnt/md0/oe_full_length/simu3/raw/output/fl_rna_simu_raw/read_tags.tsv", "r") as fs:
    line = fs.readline()
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        infos = line[0].split("_")
        if infos[1][-1] == "+":
            infos[4] = infos[4].translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdN', 'TGCAtgcaYRKMyrkmBVDHbvdhN'))[::-1]
            infos[5] = infos[5].translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdN', 'TGCAtgcaYRKMyrkmBVDHbvdhN'))[::-1]
        if line[2] == infos[4]:
            res[0] += 1
        if line[5] == infos[5]:
            res[1] += 1
        res[2] += 1
        line = fs.readline()

[114683, 114683, 114683]


res = [0, 0, 0]
with open("/mnt/md0/oe_full_length/simu3/output/fl_rna_simu/read_tags.tsv", "r") as fs:
    line = fs.readline()
    line = fs.readline()
    while line:
        line = line.strip().split("\t")
        infos = line[0].split("_")
        if infos[1][-1] == "+":
            infos[4] = infos[4].translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdN', 'TGCAtgcaYRKMyrkmBVDHbvdhN'))[::-1]
            infos[5] = infos[5].translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdN', 'TGCAtgcaYRKMyrkmBVDHbvdhN'))[::-1]
        if Levenshtein.distance(line[2], infos[4]) < 2:
            res[0] += 1
        if Levenshtein.distance(line[5], infos[5]) < 2:
            res[1] += 1
        res[2] += 1
        line = fs.readline()

[77357, 46421, 77556]




fs = open("/mnt/md0/oe_full_length/simu3/raw/fl_rna_simu_raw.fq", "r")
line = fs.readline()
while line:
    id = line.strip().split(" ")[0][1:]+"_0"
    if id not in cells_info:
        continue
    ci = cells_info[id]
    tseq = fs.readline()
    if id.split("_")[1][-1] == "-":
        bc = tseq[23:39]
        umi = tseq[39:51]
    else:
        bc = tseq[::-1][23:39].translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdN', 'TGCAtgcaYRKMyrkmBVDHbvdhN'))
        umi = tseq[::-1][39:51].translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdN', 'TGCAtgcaYRKMyrkmBVDHbvdhN'))
    if bc == ci[0]:
        res[0] += 1
    if umi == ci[1]:
        res[1] += 1
    res[2] += 1
    line = fs.readline()
    line = fs.readline()
    line = fs.readline()





