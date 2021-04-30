# Delta statistic

Codes for calculating $\delta$, a phylogenetic analog of the Shannon entropy for measuring the degreed of phylogenetic signal between a categorical trait (trait vector) and a phylogeny (metric-tree). 

## Tutorial

Make sure you have the ```ape``` package instaled before running the examples or using these scripts.

```
library(ape)
```

Upload your phylogetic tree. We exemplify with the newick format:

```
newick_tree <- "PASTE_NETWICK_TREE_HERE"
tree <- read.tree(text=newick_tree)
plot(tree)
```


It is important to guarantee that all the branches are positive as this method requires a metric-tree (i.e., branch_lengths > 0). Here, we take 1% of the 1% quantile to fill in the null branches:

```
tree$edge.length[tree$edge.length==0] <- quantile(tree$edge.length,0.1)*0.1
```

Now, we need to defined the trait vector. Confirm that the trait order follows the species order in the tree; # you can see the species order by typing: ```
tree$tip.label```.

```
trait <- c(PASTE_YOUR_TRAIT_VECTOR_HERE)
```

Now, we calculate delta:

```
deltaA <- delta(trait,tree,0.1,0.0589,10000,10,100)
```

When running the ```delta``` function you may experience this warning message: ```Warning message: In sqrt(diag(solve(h))) : NaNs produced```. Don't worry about it, it just tells me that the standard deviations of some of your rate parameters could not be calculated and these aren't used anyways.

We can also calculate p-values. Here, we shuffle the trait vector using the function delta (for 100 iterates) and create a vector of random deltas that will work as our null hypothesis. Then we compute the probability p(random_delta>deltaA) in the null distribution which returns the p-value.

```
random_delta <- rep(NA,100)
for (i in 1:100){
  rtrait <- sample(trait)
  random_delta[i] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
}
p_value <- sum(random_delta>deltaA)/length(random_delta)
boxplot(random_delta)
abline(h=deltaA,col="red")
```

* if p-value < level_of_test (generally 0.05) there is evidence of phylogenetic signal between the trait and the character
* if p-value > level_of_test there is no evidence for phylogenetic signal or the trait is saturated

etic signal

## Citation

Rui Borges, João Paulo Machado, Cidália Gomes, Ana Paula Rocha, Agostinho Antunes; Measuring phylogenetic signal between categorical traits and phylogenies, Bioinformatics, https://doi.org/10.1093/bioinformatics/bty800
