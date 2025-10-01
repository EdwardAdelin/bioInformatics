#make an app that is able to find the alphabet of a seq of text. This seq may be an ar seq or adn seq or protein.
#seq = "abbcdeefc"
#alf = [seq[0]]
for i in range(1, len(seq)):
   if seq[i] not in alf:
       alf.append(seq[i])
print(alf)


    
    