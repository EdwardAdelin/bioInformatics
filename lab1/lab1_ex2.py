#a DNA seq is given s='ACGGGCATATGCGC'. Make an app able to show precentage of the components from the alphabet of the seq S.
#in other words the input of the seq s and the output is the alphabet of the seq and the precentage of each letter
#in the alphabet found in seq s
seq = "ACGGGCATATGCGC"
alf = [seq[0]]
counts = {}
counts[seq[0]]=1
for i in range(1, len(seq)):
    if seq[i] in alf:
        counts[seq[i]]+=1
    if seq[i] not in alf:
        alf.append(seq[i])
        counts[seq[i]]=1
print(alf)
print(counts)
print("relative freq:")
total = sum(counts.values())
for key in counts:
    counts[key] = counts[key] / total * 100
print(counts)