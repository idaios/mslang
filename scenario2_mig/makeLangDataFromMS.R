library(ape)
library(phytools)
library(phangorn)
rm(list=ls())
n = 46

getJoints = function(n=NULL, intree=NULL){
    rtr = ""
    if(is.null(intree)){
        rtr = pbtree(n=n, tip.label=1:n)
    }else if(is.null(n)){
        rtr = ape::read.tree(file=intree)
    }else{
        stop("Either tree or n should be provided")
    }
    myrtr=rtr
    nl = node.depth.edgelength(rtr)
    nl[1:n] =  nl[1]
    nl = -nl + nl[1]
    myo = order(nl)[(n+1):length(nl)]
    cmdJoins=""
    for(i in myo){
        tmp = Children(rtr, node=i)
        tmp
        if( length(tmp) == 0){ next }
        flagcmd = paste("-ej ", nl[i], tmp[1], tmp[2])
        flagcmd2 = paste("-en ", nl[i], tmp[2], 0.00001)
        cmdJoins = paste(cmdJoins, flagcmd, flagcmd2)
        rtr$edge[rtr$edge[,1] == i, 1] = tmp[2]
        rtr$edge[rtr$edge[,2] == i, 2] = tmp[2]
    }
    cmdJoins
    cmdI = paste("-I",n)
    for(i in 1:n){
        cmdI = paste(cmdI, 1)
    }
    cmdI
    cmd = paste("./ms ", n, "1", cmdI, cmdJoins, "-r 100 20 -t 10 -em 0 1 40 20 -em 0 40 1 20 -T -a -seed 36585 5724 13704")
    msout = system(cmd, intern=TRUE)
    return(list(myrtr, msout, cmd))
}
####################3


##cmd="./ms 46 1 -r 10 2 -G 5 -s 512 -T -a"
simres = load("msout.RData")
simres = getJoints(n=n)##system(cmd, intern=TRUE)
msout = simres[[2]]
simtree = simres[[1]]
write(x=simres[[3]], "mscmd.txt")
##plot(simtree)

len=250

getOption = function(keyword){
    key = as.numeric(unlist(strsplit(msout[grep(keyword, msout, perl=TRUE)], split=" "))[-1])
}

getMat = function(msout){
    ind = grep(pattern="[^10]", x=msout, perl=TRUE, invert=TRUE)
    ind2 = grep(pattern="^[10]", x=msout, perl=TRUE)
    ind = intersect(ind, ind2)
    d = msout[ind]
    dmat = t(sapply(d, function(x){
        as.numeric(unlist(strsplit(x, split="")))
    }))
    row.names(dmat) = NULL
    colnames(dmat) = NULL
    dmat
}

mat = getMat(msout)

getEv = function(msout){
    probs = c(rev(0.5^(1:5)), 0.5**(1:5))
    effects = c(-5:-1, 1:5)
    names(probs) = effects
    pos = getOption("position")
    ages = getOption("ages")
    intpos = ceiling(len*pos)
    if(min(intpos) == 0){ intpos[which.min(intpos)] = 1 }
    ord = order(intpos, -ages)
    ages = ages[ord]
    intpos = intpos[ord]
    mat = getMat(msout)
    newmat = matrix(0, ncol=length(unique(intpos)), nrow=nrow(mat))
    ii = 1
    i = 1
    MAXEFFECT = 9
    while( i <= ncol(mat)){
        effect = sample(effects, 1, prob=probs)
        while( min(newmat[mat[,i] == 1, ii]) + effect < 0 || max(newmat[mat[,i] == 1, ii]) + effect > MAXEFFECT){
            effect = sample(effects, 1, prob=probs)
        }
        for(j in 1:nrow(mat)){
            newmat[j,ii] = newmat[j,ii] + effect*mat[j,i]
        }
        if(i < ncol(mat) && intpos[i] == intpos[i+1]){
            print(paste(paste(mat[,i], collapse=","), " ", effect))
            i = i+1
            next;
        }else if(i < ncol(mat) && intpos[i] != intpos[i+1]){
            print(paste(paste(mat[,i], collapse=","), " ", effect))
            print("--------------------------------")
            print(paste(newmat[,ii], collapse=","))
            print("")
            ii = ii  + 1
        }
        i = i+1
    }
    row.names(newmat) = 1:nrow(mat)
    return(newmat)
}

multMat = getEv(msout)
multMat


library(phangorn)
phy = phyDat(multMat, type="USER", levels=0:9)

getMLTree = function(phy, mat){
    treeRA = random.addition(phy)
    treeSPR = optim.parsimony(treeRA, phy)
    treeNJ = NJ(dist.gene((mat)))
    fit = pml(treeNJ, phy)
    fitGTR = optim.pml(fit, model="GTR",  optInv=T, optGamma=T, rearraangement="NNI", optEdge=TRUE, optQ=TRUE, optBf=TRUE, optRate=TRUE)
    return(fitGTR)
}

allrealTree = msout[ grep(pattern='\\(', msout)]
arealTree = allrealTree[1]
allrealTree
allrealTree[[1]]

treeCovs = as.numeric(gsub(pattern="\\[(\\d+)\\].*", replacement="\\1", x=allrealTree))
treeCovs = treeCovs/(sum(treeCovs))
if( length(treeCovs) == 0){
    treeCovs = 1
}
allrealTreeN = ape::read.tree(text=allrealTree )
realTreeN = ape::read.tree(text=arealTree[1])

write.mat = function(mat, file){
    seq = ""
    for(i in 1:nrow(mat)){
        seq = paste(seq, ">", i, sep="")
        seq = paste(seq, "\n", sep="")
        seq1 = paste(mat[i,], collapse="")
        seq = paste(seq, seq1, sep="")
        seq = paste(seq, "\n", sep="")
    }
    write(seq, file=file)
}

runraxml = function(outname, mat, model){
    write.mat(mat, outname)
    cmd = paste("raxml-ng --all --msa ", outname, " --model ", model, " --redo --bs-trees 10")
    system(cmd)
    tr = ape::read.tree(file=paste(outname, ".raxml.support", sep=""))
    tr
}

getBinary = function(multMat){
    mat = NULL
    i = 1
    jj = 1
    for(i in 1:ncol(multMat)){
        states = unique(multMat[,i])
        nstates = length(states)
        if(nstates < 2){ next }
        for(j in 1:(nstates -1 ) ){
            v = as.numeric(multMat[,i] == states[j] )
            if(is.null(mat)){
                mat = cbind(v)
            }else{
                mat = cbind(mat, v)
            }
        }
    }
    row.names(mat) = 1:nrow(multMat)
    colnames(mat) = NULL
    mat
}

binmat = getBinary(multMat)

bintree = runraxml("binmat.fa", binmat, "BIN")
multtree = runraxml("multmat.fa", multMat, "MULTI10_MK")
origbintree = runraxml("origbinmat.fa", mat, "BIN")

sizes = list(orig=ncol(mat), mult=ncol(multMat), bin=ncol(binmat))
sizes

allrealTreeN

getTreeDists = function(multTreeList, atree){
    ds = c()
    if(class(multTreeList) == "multiPhylo"){
        for(i in 1:length(multTreeList)){
            d = RF.dist(multTreeList[[i]], atree)##*weights[1]/(sum(weights))
            ds = c(ds, d)
        }
        return(ds)
    }else if(class(multTreeList) == "phylo"){
        return( RF.dist(multTreeList, atree))
    }
}

getAvDist = function(t1, t2, weights=NULL){
    if(is.null(weights)){weights = 1}
    ds = getTreeDists(t1, t2)
    avDists = sum(ds*weights/(sum(weights)))
    avDists
}


trees = list(real=allrealTreeN, guide=simtree, orig=origbintree, mult=multtree, bin=bintree)
dists = matrix(nrow=length(trees), ncol=length(trees))
colnames(dists) = names(trees)
row.names(dists) = names(trees)
for(i in 1:length(trees)){
    for(j in 1:length(trees)){
        if(i == 1 || j==1){ weights = treeCovs}
        else{ weights = NULL }
        dists[i,j] = getAvDist(trees[[i]], trees[[j]], weights)
    }
    ##dists[[i]] = dist.topo(realTreeN, trees[[i]])
}
diag(dists) = 0
dists

getTreeDists(allrealTreeN, simtree)


pdf("trees.pdf")
layout(matrix(1:6, nrow=2, byrow=TRUE))
cex=0.8
plot.phylo(simtree, main="Guide tree", cex=cex)
plot.phylo(realTreeN, main="Correct tree", cex=cex, show.node.label=F)
plot(origbintree, main="ms binary data tree", cex=cex, show.node.label=F)
plot(multtree, main="multi value tree", cex=cex)
plot(bintree, main="binarized tree", cex=cex)
dev.off()

write.table(unlist(dists), "rfdists.txt", quote=F, col.names=T, row.names=T)
write.table(unlist(sizes), "sizes.txt", quote=F, col.names=F)

save(simres, file="msout.RData")
