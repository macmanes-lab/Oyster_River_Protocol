#Response To review
## Reviewer 1

This is really interesting paper that I feel could yield good results however this manuscript feels rushed to me and I don’t think it is ready for publication yet.

Some of this is pedantic detail regarding the manuscript; the English needs to be be improved (e.g. “Troublesome, this biology is lost from all downstream analysis and inference”) as could the flow of the manuscript. There is considerable repetition (for example, at the start of the Discussion and the Conclusions) and p-values are routinely quoted to more significant figures than I believe (e.g. Section 3.2.1 p=0.0001322 for a comparison of the distributions with only 15 data-points per distribution). There are more serious shortcomings to the manuscript however that I think require considerably more attention.

> Thanks for pointing out the pedantic details related to the manuscript. I have addressed them.

I understand that the manuscript is a short presentation of a potentially useful tool but I do not feel that the authors give sufficient context in the introduction. In particular, there is no attention given to previous work on attempting to benchmark and reconcile the differences between potential de novo transcriptome assemblers, no attempt made to provide a list or or reference to the current suite of tools or any discussion of previous work highlighting how sensitive these tools are to a sensible and carefully chosen set of parameters (although the authors do correctly note that the apparent ease of use of these tools belies the complexity of the process).

> I have beefed up the introductory materials significantly.

On the plus side, the manuscript is well set out and the author has shared their code, the data they use is public and the assemblies will be made available, which is all good news.

> :)

#### Experimental design
The research question is well laid out, however the lack of summarisation of previous efforts to benchmark de novo transcriptome assembly makes it a little unclear how significant the problem they are trying to address is.

> The systematic benchmarking of assemblies, to the degree done here, is novel. Even in developing Trinity, we used Detonate (not transrate), and BUSCO did not exist at the time, nor did tools for very rapidly (and at scale) performing reciprocal blast. Therefore, I believe that these assemblies are way more benchmarked, using orthogonal metrics, that is commonplace in the assembly literature.

A good effort has been made to assess the performance of their new protocol in the context of other stand-alone assemblers, but I am concerned about several aspects of the process:

1) The main standalone tool they compare the Oyster River Protocol (ORP) to is Trinity, which itself forms one of the three assemblers that comprise the ORP. As such, it is perhaps not surprising that the ORP might outperform it.


> Thanks for the comment. I chose to compare to Trinity, mostly because this was the assembler that (arguably) most people are using. It would be the product that most people would be interested in comparing to. That the ORP outperforms Trinity is not at all a given however. It could have been possible that the merging steps introduced significant noise into the assembly, or that Trinity just outperformed the other assemblers in most ways, so contribution from other assemblers was un-helpful.

2) The author doesn’t really give any justification for why he chose these three tools for the ORP as opposed to any other of the myriad of different tools. Did he choose them because they were popular, or known to give good results, or because they have very different underlying algorithms and so might be sampling different regions of search-space?

> Excellent comment! I have included this info in the methods section. "These assemblers were chosen based on the fact that (1) use an open-science development model, whereby end-users may contribute code, (2) they are all actively maintained and are undergoing continuous development, and (3) occupy different parts of the algorithmic landscape."

> Without a doubt, other assembler could add value, but they lack these important characteristics, they were left out of the comparison. BinPacker, for instance, has produced some very nice assemblies, but it is nearly impossible to install, and the developers have been minimally responsive to my inquiries.

> Adding new assemblers, as they are developed, is definitely in the development plan. As of now, I encourage people to suggest to me different things to benchmark and add.  

3) The data preparation for the stand alone Trinity runs is not the same as for the ORP assemblies. Specifically, the data for the stand-alone Trinity runs doesn’t include the error correction step that is included in the ORP from a previously published tool (RCorrector). As a result, it is unclear whether the improvement observed for the ORP over Trinity is a result of combining the assemblies or the error correction. In particular the % read -mapping Transrate score - which is where the ORP sees the most dramatic improvement over Trinity - may well be be strongly affected by the error correction step.

> Sorry for the confusion. The reads that I used for the Trinity assembly are treated in identical fashion to the rest of the pipeline. I have clarified this in the manuscript.

4) The manuscript doesn’t give very much detail on how the BUSCO transcripts were assessed. For example, the manuscript often includes statements such as “Trinity assemblies contained on average 86% (sd = 21%) of the full-length orthologs, 154 while the ORP assembled datasets contained on average 85% (sd = 16%) of the full length transcripts.” without carefully describing what this actually means. Are these isoform matching with 100% identity over 100% of both the ISOFORM and BUSCO transcript length? I highly doubt that the matches are this perfect (from my own experience with de novo transcript assembly). We really need considerable more detail here on how these things are compared and considered matching or not.

> Here, we follow the default definitions developed as part of the BUSCO manuscript. Specifically, I have nw included in the manuscript the following language to help clarify.
"We report default BUSCO metrics as described in (Simao:2015kk). Specifically, "complete orthologs", are defined as query transcripts that are within 2 standard deviations of the length of the BUSCO group mean, while contigs falling short of this metric are listed as "fragmented"."

> BUSCO definitions are not that tight with regards to ISOFORMS, and, on an aside, proper isoform reconstruction using Illumina data is not something to be trusted, as you are aware. I agree there is significant complexity here, but in the absence of real confidence in isoform reconstruction, and a well benchmarked tool (BUSCO), I think that reporting default values is helpful to readers.


5) The inclusion of multiple kmer lengths is interesting, but no real effort is made to determine whether the differences between the individual assemblies is due to the different kmer lengths or the different algorithms of the tools. For example, how does an ORP assembly compare with three Trinity assemblies, with kmers of 25,55 & 75, combined in the same fashion?

> This is impossible, as the longest kmer that Trinity can use is k=32. I have more fully fleshed out the contributions of k=55 and k=75 within SPAdes. See figures 3, 4, and 5, as well as Table 2. This (especially the MDS plot and table 2) demonstrate value in running SPAdes at the 2 different kmer lengths.    

Overall, however, the manuscript is lacking the depth and detail I think it needs to really fill the knowledge gap and without this details reproducing these results would be difficult despite the easy and availability of the software install and the dataset.

> I have added significant detail in all phases of the manuscript which should help in the interpretation of the results. Regarding reproducibility, all the data are publicly available, and the code I used is freely available on Github.

#### Validity of the findings
In its current form, the manuscript has not convinced me that the ORP represents any significant improvement over Trinity, the existing stand-alone de novo assembly tool tested. As such, I think that the conclusion the the ORP represents an improvement over existing tools is not yet sound.

> I hope that the revised manuscript, with the additional benchmarking, has convinced you.

#### Comments for the Author
This work details a new multi-assembler approach to improving de novo transcriptome assembly. This is a really good idea that is worth exploring and I really want to see a strong detailed exploration of the concept and how much we might gain from combining de novo transcriptome tools, but I’m afraid that these have not been adequately illuminated here. A little more time, exploration of the tools and more careful analysis and I think it would make a good paper. I'd be happy to review future revisions of the paper and would like to see ORP broadened to include more assembly tools.

> There is a balance between adding more tools and the time it takes to run the ORP. I have tried to justify the choice of assemblers better (Lines 105-107). If there are specific assemblers that meet those criteria that you'd like me to try, please suggest them specifically!

The instruction for installing and running the software (http://oyster-river-protocol.readthedocs.io/en/latest/aws_setup.html) work and I could set up the tool on my local box and verify that it runs. Unfortunately I don't have access to sudo on my clusters large memory machines, and we have no way of charging AWS tie to grants, so I couldn't test it on live data. Perhaps, given the complex web of dependencies involved, it might be an idea to distribute this as a conda package or a docker instance?

> I will consider making a bioconda/docker image, but I do not have the bandwidth at the moment (I've never developed a bioconda package). Open to collaborations on this! I will say, however, that the install (once you have Linuxbrew) is easy (`make` in install dir), and does not require `sudo` privileges. I have stated this in the manuscript (line 85), and have reiterated in the online instructions.

## Reviewer 2
#### Basic reporting
The manuscript is well written and well organized.

The transcriptome assemblies resulting from each of the assemblies are made available as supplementary data, but there are many other supplementary data files and figures that would be useful to include that were not available, such as the TransRate reports for each of the assemblies, and ideally a summary data table that includes all TransRate scores and score components (ie. the individual metrics computed by TransRate that are leveraged in computing the final overall assembly score). Busco scores for all assemblies should be included as well.

> I have compiled these results into a csv file, available at https://github.com/macmanes-lab/Oyster_River_Protocol/blob/master/manuscript/orp.csv

Although multiple transcriptome assemblers were leveraged, the figures focused almost entirely on comparisons between the single assembler Trinity and the author's results from ORP. It would be useful to include comparisons to the other assemblers as well in a similar format, and provided as supplementary materials if not deemed worthy as main figures.

> Excellent idea. I have made available all the data available in a spreadsheet, so any number of comparisons (too many to reasonably include) are possible. Of course I compare the ORP to Trinity, as it's widely regarded as the best, and is the most used of the group. In many places in the manuscript, while still focusing on trinity, I now report information related to the performance of the other assemblers. In no case do these results contradict the main findings, though in some cases the differences between the ORP assembly and SPAdes/Shannon are more or less than the differences between ORP and Trinity. Does this make sense?

#### Experimental design
Comparisons among the final assemblies involve the TransRate final assembly scores, percent of reads represented by each assembly, and core ortholog representation as per Busco. These are all useful metrics, but it would be of interest to include additional metrics that reflect rates of chimerism.

> Thanks Brian, chimerism is indeed a big problem, and a hard one to directly detect, particularly when lacking a reference. TransRate does provide a metric related to chimerism, and I now report that metric (Line 202-206). Rates of chimerism range between 10% and 12%, with Shannon being the lowest (ORP is the 2nd lowest at 10.4%). So, what I can say is that the ORP does not improve or worsen chimerism.  

It would benefit the manuscript to include an orthogonal measure of assembly quality such as provided by the DETONATE software. Ideally, this would show that ORP compiled assemblies have improved DETONATE scores as compared to the assemblies provided as input to ORP.

> Thanks Brian for pushing me to do this. I have added the detonate scores to the analyses (See figure 1c). Without exception, the scores are improved in the ORP assembly over the Trinity assemblies.

While ORP attempts to stringently cluster transcripts into isoform-level groupings as a way of minimizing the loss of relevant splicing isoforms, there's no evidence shown of the impact of this method on splice isoform representation in the final ORP assembly. This deserves much more attention in the manuscript. Also, in the case of isoforms that are fused as chimeras, selecting a single best representation from a cluster containing chimeras and non-chimeras would expect to contribute to loss of important transcripts. For example, given a cluster {AB, A, B} where AB is a chimera, selecting a single representative transcript with the best score would presumably yield A or B and hence exclude an important transcript in the final output. More attention to resolving clusters of transcripts containing chimeras is warranted.

> I agree that this ideally deserves more attention, but I just don't see a way to do that, given that all but 2 of these assemblies lack high-quality references. In addition, it is my long-standing belief (I think we agree), that isoform reconstruction is inherently flawed in short-read assemblies. We just have such limited ability to accurately reconstruct many/most splice isoforms, and in most cases can't even know which ones are right and which ones are incorrectly assembled. I have added some discussion to this effect in the paper (lines 217-224), to warn readers. In short, the way contig scores are calculated, coupled with the improvement of the scores related to assembly structure seem to suggest that this type of loss is not common.

Since ORP is centered around combining results from multiple assemblers, it begs the question of what assemblers would be generally best suited to use with ORP. The author leverages Trinity, Shannon, and Spades. Trinity and Shannon have been previously demonstrated to be highly effective transcriptome assemblers, but to my knowledge, Spades has seen comparably little use for transcriptome assembly; instead Spades has been more targeted towards genome assembly. There does appear to be an rnaSPAdes method included with Spades that's more specifically targeted to transcriptome assembly, but I could find little information about it being used in practice nor how effective it is for resolving spliced isoforms, which is one of the primary goals of transcriptome assembly and one that most differentiates it from other assembly challenges. One of the goals of ORP is to leverage multiple kmer lengths as part of the assembly, and Spades is run twice using different kmer values as one way to accomplish this. Using multiple kmer values to generate an assembly, leveraging the benefits of small and large kmers, is a feature of existing transcriptome assemblers Oases and IDBA-tran, neither of which were used by ORP nor cited. It isn't clear as to why the author would rely so heavily on Spades when other highly effective transcriptome assemblers are readily available and already implement ideas that are central to ORP's goals. ORP would ideally be shown to benefit from using the very best assemblers as well as to be resilient towards including results from assemblers that are sub par.

> I have added some descriptions of the rationale I used for selecting assemblers. In short, they had to be active developed and open sourced. I benchmarked a larger number of assemblers, and picked the ones that seemed to come from different parts of the "algorithmic landscape". I also picked ones that I could reliably install and run (BinPacker was excluded for this reason, even though when it worked, it performed very well). I describe this decision making process in lines 112-115.

> I have beefed up the into, to include references to some more of the classical assemblers, including Oases and IDBA-trans, which indeed both made important earlier contributions to the field.


#### Validity of the findings
The primary finding is that the author's ORP yields assemblies that are generally higher in quality than the Trinity assembler alone, primarily leveraging TransRate and BUSCO for assembly quality assessment. In terms of the bulk statistics, this is a valid conclusion, but it does not specifically address reducing chimeric contigs and ensuring representation of alternatively spliced isoforms, both very important aspects of transcriptome assembly and also often negligible with respect to the bulk statistics focused upon here.

> I think I covered this above, but just to recap. I've included metrics related to chimerism (they are stable across all assemblers) and provide some enhanced discussion about splice isoforms in the discussion. Here, I argue that the (indirect) measures suggest this is not a huge problem, but given the lack of reference genomes and generally poor-performance of short-reads in this specific area, its hard to really know for sure.   


Overall, the results section is overly directed towards comparisons to Trinity when there were several other assembly methods being used by ORP. For each assembly quality metric, it would be useful to see how ORP compares to each of the individual assemblies, including Trinity. If there exist cases where any individual assembly used as input to ORP ends up outperforming ORP, it would call into question the value of applying ORP as a general protocol.

> I've added language throughout where I note the performance relative to the other assemblies. Other assemblers underperform the ORP, though in a very small number of cases, they outperform Trinity. I still think that the major comparison needs to be to Trinity as this is the most commonly used assembler, and the thing most people will be familiar with. Does this make sense?


#### Additional comments:

As a protocol, ORP should apply equally well to de novo transcriptome assemblies as to genome-guided assemblies. Using reference model organism data sets, one could explore applying ORP to transcript sets resulting from cufflinks, stringtie, or other genome-guided methods such that ORP will yield the best quality transcript set.

> I agree that leveraging genome-guided approaches could be beneficial. This is something I'd like to explore in the future, but seems outside of the scope of the current work.

Most data sets in the manuscript appear to lightly benefit from ORP, whereas there are a few assemblies that greatly benefit from ORP, especially in the context of percent reads mapping. It would be useful to comment on this and whether the suboptimal Trinity assemblies in these cases were reflective of some aspect of the data, and whether the alternative assemblers provided much better assemblies than Trinity in these cases. This would be much easier to assess if the results from all methods were easily accessible in tables and figures in supplementary materials.

> True, most assemblies are improved by ~2-10% over single assembler solutions, across several metrics, so the improvement is marginal. I think given the modest increase in computational requirements (and 0 increase in hands on time), this improvement will be valuable to researchers wi=anting to squeeze any last drop of utility out of their assemblies. About the specific contrast to Trinity, sometimes yes the Trinity mapping is the worst, and sometimes not. I've added language in multiple places in the revised manuscript to illustrate the relationships between the ORP and other assemblers.  

The composition of the final ORP assemblies appears to mostly derive from Spades assembled transcripts. I imagine this could happen for several reasons. For example, (a) there could be clusters of transcripts where Spades yields the best representative transcript with the highest score, (b) there are clusters where there are tied scores and Spades is chosen preferentially, (c) Spades generates many transcripts that are unique to Spades and these are automatically propagated to ORP. It would be useful to know how the composition of the final ORP assembly reflects which transcripts were selected as best within clusters vs. alternative justifications for their inclusion in the final assembly.

> I think this was largely an artifact of how I described the composition in the last version of the manuscript. Basically, each assembly contributes about 25% of the final assembly. So technically SPAdes does contribute more to the final assembly that do the other 2, but I think this is largely a function of the fact that we use 2 different kmer lengths.

> About the alternative explanation B. I looked a few hundred scored clusters, and no ties were observed (at the high-scoring end). The scores are reported to 5 significant digits, so a lot of room to be different. In the case of ties, the 1st contig is selected in the list or transcripts. This list comes from OrthoFinder, and is populated in random order (per David Emms, the Orthofinder developer). So, I don't think B is likely. Alternative C: I've made a new Table 2 which reports the number of unique transcripts, which is a very small number relative to the total number of transcripts. I also wondered if it might be due to the fact that there are just more SPAdes transcripts to chose from (higher contig number). In some cases that number of SPAdes transcripts exceeds other assembles (this is try especially for SPAdes kmer=55), this is not the case (https://github.com/macmanes-lab/Oyster_River_Protocol/blob/master/manuscript/contrib_by_numcontigs.eps). There is no relationship between the number of contigs in the primary assembly and the percent contribution.

## Reviewer 3
#### Basic reporting


In this paper, Dr. MacManes describes a new computational protocol for de novo transcriptome protocol. Briefly, the protocol first runs three different assemblers (Trinity, Shannon, and SPAdes) with multiple parameters, combines the resulting assemblies and extracts ortholog groups using a custom modification of OrthoFinder called OrthoFuse, and finishes using a custom version of Transrate to select the highest contig score from each ortholog group. Dr. MacManes demonstrates that there are very few downsides to using this protocol over e.g. a standard Trinity protocol: many metrics are either the same or improved.

This is an important paper for a wide range of biologists, because (as is correctly noted in the paper) de novo transcriptome assembly is very widely used in research into the biology of many (most!) different eukaryotic organisms. The completeness of reference transcriptomes is important for all the downstream analyses, including annotation, recovery of homologs, and inference of differentially expressed genes.

Along these lines, there are a number of conclusions in the paper that are quite interesting and important:
* mapping rates and transrate assembly scores are significantly improved in the ORP assemblies.
* Busco completeness scores are generally the same or improved in the ORP assemblies.
* the different assemblers appear to be better at recovering transcripts from different abundance fractions of the reads (which I find really interesting and unexpected, at least in the degree to which it seems to occur). This observation will likely lead to methodological advances in standalone assemblers, which is excellent.

The work in general seems to be of high quality and many of the conclusions are solidly grounded in the results and discussion.

However, there are a number of puzzling choices and omissions in the paper as it stands. I detail these below. There are also many minor typos and grammatical errors that I will communicate directly to the author.

#### Questions and major comments:

* line 52 - I would imagine that other kinds of biases might result from the approach chosen in the paper; any thoughts as to downsides of the ORP approach here? (See comment on evaluating against a quality ref transcriptome, below.) e.g. in genome assemblies, maximizing N50 leads to misassemblies; might that happen here?

> reviewer 2 brings up the same points, specifically as it relates to isoforms and chimerism. While I don't have a great answer to either issue, mainly because of the lack of high-quality reference genomes, I can use the TransRate segmentation metric to estimate chimerism. I trust short-read isoform assembly so little that I don't have much insight here. I added in some discussion and indirect reasoning about why I don't think the ORP is doing worse than other things, but it's probably not doing much better also.

* unweighted assembly content is not systematically evaluated in this paper - most of the evaluations such as mapping rate and transrate score are biased by underlying abundance distribution of the transcripts, i.e. more abundant transcripts will produce more reads. I would suggest doing some kind of unweighted assembly comparison (one with higher resolution than BUSCO; next point) to see if more overall content is produced. This could be what Figure 4 intends to show, at least partly, but see below.

> I now use Shmlast analysis against Swissprot for a more higher resolution comparison. See Figure 3 and Table 2 for a summary. In all cases, the scores are better for the ORP reconstruction. Importantly, this analysis helped me identify a place where further optimization can be had. Specifically, the ORP seems to leave out some real transcripts that 1 of the individual assemblers reconstructed. I have a student working on tracking these transcripts down, and this will be included in a future ORP release. 

* it seems that no evaluation is done against a quality reference transcriptome, e.g. mouse RNAseq -> mouse. This is a major omission; the major reference-based approach used, BUSCO, only looks at a small subset of the transcripts, and does not evaluate transcript content systematically. It would be valuable to use a reference-based measure of transcriptome quality here, such as the Trinity paper used.

> Thanks for the important suggestion. I have attempted to address this by using the `Shmlast
` package to find CRBH's in the Swissprot database, which I believe is better that looking in individual (e.g. Mus) reference transcriptomes, in that the same reference database can be used for all assemblies, thereby allowing for inter-dataset comparison. What this analysis shows, is that in all cases, I recover more CRHBs using the ORP that when using any other assembly, which speaks to the general merits of the approach.

* The choice of taxonomic descriptions in table 1 is bewildering - I don't see any reason (biological or other) for the descriptions chosen. Maybe stick with 'animal' and 'plant'?

> I reworked the figure to include animal, plant, and protozoa. I agree this makes the presentation much cleaner.

* More, section 4.2 incompletely addresses an interesting question - there are many features of genomes/transcriptomes that are more likely to challenge transcriptome assembly than taxonomic membership, e.g. tissue sampled, tissue complexity in terms of number of cell types, presence of RNA editing, recent genome duplications, polymorphism rate, etc. It's probably beyond the scope of this paper to explore these issues, but I'm not sure how it merits the incomplete discussion it's given in the current paper. Perhaps a single sentence "ORP performed equally well across all of the organisms chosen at random" would be better, or something like that.

* The paper states that no value is seen in including reads beyond 30m, which is intuitive and matches our experience. However, most transcriptome assembly projects sequence multiple tissues and then combine them in a de novo assembly. Is the recommendation for the ORP to subsample each tissue data set to 30m and then combine them? Or what? Some comment (along with a specific reference to a paper with an overall suggested workflow for making use of the assembly) would be welcome.

 > My general recommendation for people strolling on to my office is to subsample each tissue set to ~30M pairs, assemble, and merge. If isoforms are of zero interest (or you have zero confidence in isoform assembly from short read data, as I do), then the recommendation is to combine reads from all tissues together then generate one 30M read dataset. The assembly of this dataset is likely to contain some isoforms not seen in "real life" (e.g., imagine chimerism of liver-specific isoform with brain-specific isoform), but given isoform reconstruction is so bad anyway, the risk of these chimeric events might not be so impactful.

 > I'd be *really* happy to provide some comment in the manuscript, but I'm afraid that the recommendation is based on a lot of anecdote, and very little data. Further, I don't really know of a paper where this question has been fully resolved. Do you know of one?

* Along those lines, the extra computational requirements of the ORP seem likely to be significant when combining multiple tissue RNAseq data sets. This is one of the reasons that Trinity recommends using in silico normalization, which was explicitly *not* used in this paper. The discussion should be expanded to include this.

> I generally do not recommend diginorm, as it (in my hands) has resulted in slightly less contiguous assemblies. I can certainly add this  

* I'm puzzled by the division between Results and Discussion in this paper. I'd suggest either combining them or having the primary results be separate from the discussion of what the results mean; typically the division I use is "here are the facts I see" (for results), and "here is what I think they mean".

* I had trouble understanding the section starting at line 188, and figure 4, at least at first. If I understand correctly, this paragraph is discussing both swissprot match percentages AND transcript "origin" percentages. The latter can be studied because orthofuse doesn't actually merge transcripts, it picks a "best" transcript from each orthologous grouping - is that right? The former (swissprot analysis) is of objective interest because it speaks to what genes are being recovered by which assembler. Or... maybe that's not what the paragraph is saying. It would be nice if this could be clearly split into discussions of homology/orthology/recover of known genes, and recovery of transcripts (known or unknown). As it is I think it's muddled in its message. Figure 4 seems clear enough but I'm not sure if it's talking about transcripts (my guess) or swissprot matches. (In retrospect, it seems like the beginning of section 3.2.3 is talking about transcript origin wrt the final assembly, while the bottom of that section is talking about swissprot).

#### Minor comments:

* surely there is a review citation of some sort for the statement that "transcriptome sequencing has been influential"? p2 line 24
* Throughout the paper, clickable links in the PDF are used (The code ... is available _here_.) these should be replaced with URLs, perhaps in footnotes or references.

> These have all be replaced.
