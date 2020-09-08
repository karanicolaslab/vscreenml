import re
import sys

file_name = sys.argv[1]
f = open(file_name, "r")
content = f.read()
f.close()
data = {}

for block in content.split("\n\n"):
    block = block.lstrip()
    block = block.rstrip()
    # print block

    if re.search("Atom-type pair counts within 2.5 angstroms:", block):
        for l in block.split("\n"):
            matchObj = re.search(r"Atom|.*----", l)
            if not matchObj:
                k, v, p = l.split("|")
                k = k.strip()
                v = v.strip()
                p = p.strip()
                # print (k,v,p)
                k_v = str(k) + "_" + "2.5A" + "_" + str(v)
                # print k_v
                # print float(p)
                data[k_v] = float(p)
            # print data[k_v]

    if re.search("Atom-type pair counts within 4.0 angstroms:", block):
        for l in block.split("\n"):
            matchObj = re.search(r"Atom|.*----", l)
            if not matchObj:
                k, v, p = l.split("|")
                k = k.strip()
                v = v.strip()
                p = p.strip()
                # print (k,v,p)
                k_v = str(k) + "_" + "4.0A" + "_" + str(v)
                # print k_v
                # print float(p)
                data[k_v] = float(p)
                # print data[k_v]

    if re.search("Summed electrostatic energy by atom-type pair", block):
        TotalElec = 0.0
        for l in block.split("\n"):
            matchObj = re.search(r"Summed|Atom|.*----", l)
            if not matchObj:
                k, v, p = l.split("|")
                k = k.strip()
                v = v.strip()
                p = p.strip()
                # print (k,v,p)
                k_v = str(k) + "_" + "Elec" + "_" + str(v)
                # print k_v
                # print float(p)
                TotalElec += float(p)
        data["TotalElec"] = TotalElec
        # print data["TotalElec"]

    if re.search("Number of rotatable bonds in the ligand:", block):
        for l in block.split("\n"):

            v, p = l.split(":")
            v = v.strip()
            p = p.strip()
            # print (v,p)
            k_v = "RotBond"
            # print k_v
            # print float(p)

            data[k_v] = float(p)

    if re.search("Active-site flexibility:", block):
        for l in block.split("\n"):
            matchObj = re.search(r"Active-site|Sidechain/Backbone|.*----", l)
            if not matchObj:
                k, v, p = l.split("|")
                k = k.strip()
                v = v.strip()
                p = p.strip()
                # print (k,v,p)
                k_v = str(k) + "_" + "flex" + "_" + str(v)
                # print k_v
                # print int(p)
                data[k_v] = float(p)
                # print data[k_v]

    if re.search("Hydrogen bonds:", block):
        TotalHbond = 0.0
        for l in block.split("\n"):
            matchObj = re.search(r"Hydrogen bonds|Location of Donor|.*----", l)
            if not matchObj:
                r, k, v, p = l.split("|")
                r = r.strip()
                k = k.strip()
                v = v.strip()
                p = p.strip()
                # print (r,k,v,p)
                k_v = str(r) + "_" + str(k) + "_" + str(v) + "Hbond"
                # print k_v
                # print float(p)
                TotalHbond += float(p)
        data["TotalHbond"] = TotalHbond
        # print data["TotalHbond"]

    if re.search("Hydrophobic contacts", block):
        TotalHphob = 0.0
        for l in block.split("\n"):
            matchObj = re.search(r"Hydrophobic contacts|Sidechain/Backbone|.*----", l)
            if not matchObj:
                k, v, p = l.split("|")
                k = k.strip()
                v = v.strip()
                p = p.strip()
                # print (k,v,p)
                k_v = str(k) + "_" + "Hphob" + "_" + str(v)
                # print k_v
                # print float(p)
                TotalHphob += float(p)
        data["Hphob"] = TotalHphob
        # print data["Hphob"]

    if re.search("pi-pi stacking interactions:", block):
        TotalPipi = 0.0
        for l in block.split("\n"):
            matchObj = re.search(r"pi-pi|Secondary Structure|.*----", l)
            if not matchObj:
                v, p = l.split("|")
                v = v.strip()
                p = p.strip()
                # print (v,p)
                k_v = str(v) + "_" + "Pipi"
                # print k_v
                # print float(p)
                TotalPipi += float(p)
        data["Pipi"] = TotalPipi
        # print data["Pipi"]

    if re.search("T-stacking (face-to-edge) interactions:", block):
        TotalTstac = 0.0
        for l in block.split("\n"):
            matchObj = re.search(r"T-stacking|Secondary Structure|.*----", l)
            if not matchObj:
                v, p = l.split("|")
                v = v.strip()
                p = p.strip()
                # print (v,p)
                k_v = str(v) + "_" + "Tstac"
                # print k_v
                # print float(p)
                TotalTstac += float(p)
        data["Tstac"] = TotalTstac
        # print data["Tstac"]

    if re.search("Cation-pi interactions:", block):
        TotalCapi = 0.0
        for l in block.split("\n"):
            matchObj = re.search(r"Cation-pi|Secondary Structure|.*----", l)
            if not matchObj:
                k, v, p = l.split("|")
                k = k.strip()
                v = v.strip()
                p = p.strip()
                # print (k,v,p)
                k_v = str(k) + "_" + "Capi" + "_" + str(v)
                # print k_v
                # print float(p)
                TotalCapi += float(p)
        data["Capi"] = TotalCapi
        # print data["Capi"]

    if re.search("Salt Bridges:", block):
        TotalSalt = 0.0
        for l in block.split("\n"):
            matchObj = re.search(r"Salt Bridges:|Secondary Structure|.*----", l)
            if not matchObj:
                v, p = l.split("|")
                v = v.strip()
                p = p.strip()
                # print (v,p)
                k_v = str(v) + "_" + "Salt"
                # print k_v
                # print float(p)
                TotalSalt += float(p)
        data["Salt"] = TotalSalt
    # print data["Salt"]
# print "SIDE_flex_ALPHA","SIDE_flex_BETA","SIDE_flex_OTHER","BACK_flex_ALPHA","BACK_flex_BETA","BACK_flex_OTHER","TotalElec","TotalHbond","Hphobe_contact","Pi_Pi","T-stack","Cation-pi","Salt_Bridge","RotBonds"
print(
    data.get("SIDECHAIN_flex_ALPHA", 0.0),
    data.get("SIDECHAIN_flex_BETA", 0.0),
    data.get("SIDECHAIN_flex_OTHER", 0.0),
    data.get("BACKBONE_flex_ALPHA", 0.0),
    data.get("BACKBONE_flex_BETA", 0.0),
    data.get("BACKBONE_flex_OTHER", 0.0),
    data.get("TotalElec", 0.0),
    data.get("TotalHbond", 0.0),
    data.get("Hphob", 0.0),
    data.get("Pipi", 0.0),
    data.get("Tstac", 0.0),
    data.get("Capi", 0.0),
    data.get("Salt", 0.0),
    data.get("RotBond", 0.0),
)
