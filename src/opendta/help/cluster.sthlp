{smcl}
{title:Title}

{p 4 4 2}
{bf:cluster} {hline 2} Cluster analysis

{title:Syntax}

{p 8 16 2}
{cmd:cluster} {cmd:singlelinkage} | {cmd:averagelinkage} | {cmd:completelinkage} | {cmd:waveragelinkage} | {cmd:centroidlinkage} | {cmd:medianlinkage} | {cmd:wardslinkage} varlist [if] [in] [{cmd:,} {cmd:name(}clname{cmd:)} {cmd:measure(}measure{cmd:)}]

{p 8 16 2}
{cmd:cluster} {cmd:kmeans} | {cmd:kmedians} varlist [if] [in] {cmd:,} {cmd:k(}#{cmd:)} [{cmd:name(}clname{cmd:)} {cmd:measure(}measure{cmd:)} {cmd:start(}firstk|lastk|random(#)|everykth{cmd:)} {cmd:iterate(}#{cmd:)} {cmd:generate(}newvar{cmd:)} {cmd:keepcenters}]

{p 8 16 2}
{cmd:cluster} {cmd:generate} newvar {cmd:=} {cmd:groups(}#{cmd:)} | {cmd:cut(}#{cmd:)} [{cmd:,} {cmd:name(}clname{cmd:)} {cmd:ties(}error|skip{cmd:)}]   {cmd:cluster} {cmd:dendrogram} [clname] [{cmd:,} {cmd:cutnumber(}#{cmd:)} {cmd:cutvalue(}#{cmd:)} {cmd:labels(}varname{cmd:)} {cmd:showcount} {cmd:horizontal}]

{p 8 16 2}
{cmd:cluster} {cmd:stop} [clname] {cmd:,} {cmd:rule(calinski}|{cmd:duda)} {cmd:groups(}numlist{cmd:)}   {cmd:cluster} {cmd:list} | {cmd:dir} | {cmd:drop} | {cmd:rename} ...

{p 8 16 2}
{cmd:matrix} {cmd:dissimilarity} matname {cmd:=} [varlist] [if] [in] [{cmd:,} {cmd:measure} {cmd:observations} {cmd:variables}]

{title:Description}

{pstd}
Hierarchical clustering uses the Lance-Williams recurrence on the dissimilarities; ties between equal dissimilarities are broken by the smallest pair (VERIFY). K-means and k-medians use Lloyd's algorithm from the chosen starting centers. The measures are L2, L2squared, L1, Linf, Lp, Canberra, correlation and angular, and for binary data (nonzero is 1) the similarities matching, Jaccard, Russell, Hamann, Dice, antiDice, Sneath, Rogers, Ochiai, Yule, Anderberg, Kulczynski, Pearson and Gower2 (clustered as 1 - similarity). {cmd:clustermat} {it:linkage} matname [{cmd:,} {cmd:shape(}full|lower|llower|upper|uupper{cmd:)} {cmd:add} {cmd:clear} {cmd:labelvar(}newvar{cmd:)} {cmd:name()} {cmd:generate()}] clusters a dissimilarity matrix. {cmd:keepcenters} appends the final group centers to the data. Cluster information is stored in the dataset characteristics, so it is saved with the data.

{pstd}
Not implemented yet: the {cmd:dissim()} option of {cmd:matrix dissimilarity}.
