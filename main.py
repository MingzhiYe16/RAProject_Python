import csv
import random

# Global parameter
PoolSize=30000
SampleSize=50000
errorProb=0.01
AmiNoAcids=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

# Generate sample pool
def getCDR3():
    res=random.choices(AmiNoAcids,k=15)
    return ''.join(res)
Original_Alpha_Pool=[getCDR3() for i in range(PoolSize)]
Original_Beta_Pool=[getCDR3() for i in range(PoolSize)]
SamplePoolAlpha=random.choices(Original_Alpha_Pool,k=PoolSize)
SamplePoolBeta=random.choices(Original_Beta_Pool,k=PoolSize)
SamplePool=list(zip(SamplePoolAlpha,SamplePoolBeta))

IndexnAnBMap=[]
with open('IndexnAnBMap.csv', 'r') as csvfile:
    reader = csv.reader(csvfile, skipinitialspace=True)
    for row in reader:
        IndexnAnBMap.append(row)
IndexnAnBMap=[[int(x[0]),int(x[1])] for x in IndexnAnBMap]

def getSample():
    NumOfChains=random.choices(IndexnAnBMap,k=SampleSize)

    def getChainTypes(i):
        res=[]
        for j in range(NumOfChains[i][0]):
            res.append('TRA')
        for j in range(NumOfChains[i][1]):
            res.append('TRB')
        return res

    chains=[getChainTypes(i) for i in range(PoolSize)]
    chains = [item for sublist in chains for item in sublist]

    def getBarcode(i):
        res=[]
        for j in range(sum(NumOfChains[i])):
            res.append('cell_'+str(i+1))

        return res

    barcodes=[getBarcode(i) for i in range(PoolSize)]
    barcodes = [item for sublist in barcodes for item in sublist]

    def getCDR3sForCertainClonotype(i):
        result=[[],[]]
        numberA,numberB=NumOfChains[i][0],NumOfChains[i][1]
        indexInPool = random.sample(range(PoolSize), 3)
        usedA,usedB=set(),set()
        if numberA >= 1:
            if random.uniform(0,1) < errorProb:
                result[0].append(SamplePool[indexInPool[0]][0])
                usedA.add(indexInPool[0])
            else:
                result[0].append(SamplePool[indexInPool[1]][0])
                usedA.add(indexInPool[1])
        numberA-=1

        if numberB >= 1:
            if random.uniform(0,1) < errorProb:
                result[1].append(SamplePool[indexInPool[0]][1])
                usedB.add(indexInPool[0])
            else:
                result[1].append(SamplePool[indexInPool[2]][1])
                usedB.add(indexInPool[2])
        numberB-=1

        while numberA > 0:
            indexInPool=random.randint(0, PoolSize-1)
            while indexInPool in usedA:
                indexInPool = random.randint(0, PoolSize - 1)
            usedA.add(indexInPool)
            result[0].append(SamplePool[indexInPool][0])
            numberA-=1

        while numberB > 0:
            indexInPool=random.randint(0, PoolSize-1)
            while indexInPool in usedB:
                indexInPool = random.randint(0, PoolSize - 1)
            usedB.add(indexInPool)
            result[1].append(SamplePool[indexInPool][1])
            numberB-=1

        return result[0]+result[1]

    CDR3=[getCDR3sForCertainClonotype(i) for i in range(PoolSize)]
    CDR3 = [item for sublist in CDR3 for item in sublist]

    sample=list(zip(barcodes,chains,CDR3))
    sample=[('barcode','chain','cdr3')]+sample
    return sample
sample=getSample()
with open("sample.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(sample)