import rnarryext
import random

d = {}
for i in range(100):
    for c in 'ATGC':
        r = list(c * 172 + 'ATGC-') * 17
        random.shuffle(r)
        d[i, c] = ''.join(r)

g = rnarryext.RandomReadCounter(d)
#print g.simulate('ATTGGCAGAGCGTGCTTAGC', 1000000)
#print g.simulate('ATGCAAGCTAGTGAAGTCAG', 100)
r = g.simulate('TTTTGAGGCTTTCTCCCAACGCACAGACTTGTGTAATTCTAACACTAATCCTGTGAAGGGTTGTGGTTGACAGCTGGAGCCTGGGTGACATTCTACA', 1000)
print r
print map(rnarryext.shannon_entropy, r)

