{smcl}
{title:Title}

{p 4 4 2}
{bf:ca}, {bf:mca}, {bf:mds}, {bf:mdsmat}, {bf:procrustes}, {bf:biplot} {hline 2} Correspondence analysis, multidimensional scaling, Procrustes analysis and biplots

{title:Syntax}

{p 8 16 2}
{cmd:ca} rowvar colvar [if] [in] [fweight] [{cmd:,} {cmd:dimensions(}#{cmd:)} {cmd:normalize(}symmetric|standard|row|column|principal{cmd:)}]

{p 8 16 2}
{cmd:camat} matname [{cmd:,} {cmd:rowname(}string{cmd:)} {cmd:colname(}string{cmd:)} {cmd:dimensions(}#{cmd:)} {cmd:normalize()}]

{p 8 16 2}
{cmd:mca} varlist [if] [in] [{cmd:,} {cmd:dimensions(}#{cmd:)} {cmd:normalize(}adjusted|principal{cmd:)}]

{p 8 16 2}
{cmd:mds} varlist [if] [in] [{cmd:,} {cmd:id(}varname{cmd:)} {cmd:dimension(}#{cmd:)} {cmd:measure(}measure{cmd:)} {cmd:method(}classical|modern{cmd:)}]   {cmd:mdsmat} matname [{cmd:,} {cmd:dimension(}#{cmd:)} {cmd:method()}]   {cmd:mdslong} depvar [if] [in] {cmd:,} {cmd:id(}var1 var2{cmd:)} [{cmd:dimension(}#{cmd:)} {cmd:method()} {cmd:force}]

{p 8 16 2}
{cmd:procrustes} ({it:target varlist}) ({it:source varlist}) [if] [in] [{cmd:,} {cmd:transform(}orthogonal|unrestricted{cmd:)} {cmd:norho}]

{p 8 16 2}
{cmd:biplot} varlist [if] [in] [{cmd:,} {cmd:std} {cmd:dim(}# #{cmd:)} {cmd:rowlabel(}varname{cmd:)} {cmd:norow} {cmd:nocolumn}]

{title:Description}

{pstd}
{cmd:ca} analyzes a two-way table by singular value decomposition of the standardized residuals. {cmd:mca} analyzes the indicator matrix of several variables, with Benzecri-adjusted inertias by default (VERIFY). {cmd:mds} and {cmd:mdsmat} scale a dissimilarity matrix: classical scaling by the eigenvalues of the double-centered squared dissimilarities, modern scaling by stress majorization from the classical solution. {cmd:camat} analyzes a matrix of counts. {cmd:mdslong} scales dissimilarities stored one pair of objects per observation. {cmd:procrustes} finds the orthogonal rotation, dilation and translation (or, with {cmd:transform(unrestricted)}, the least-squares linear transformation) that bring the source variables closest to the target. {cmd:biplot} draws the observations and the variables on the first two singular dimensions of the centered data.

{pstd}
Not implemented yet: supplementary rows, columns and variables, {cmd:mca} with {cmd:method(indicator)} and {cmd:method(johnson)}, {cmd:mds} with {cmd:method(nonmetric)} and {cmd:loss()}, {cmd:mdslong} with {cmd:s2d()}, {cmd:procrustes} with {cmd:transform(oblique)}.

{title:Postestimation}

{phang2}{cmd:estat} {cmd:table} | {cmd:coordinates} (after {cmd:ca} and {cmd:mca}); {cmd:estat} {cmd:config} (after {cmd:mds}){p_end}
{phang2}plots: {cmd:cabiplot} and {cmd:caprojection} (after {cmd:ca} and {cmd:camat}; options {cmd:dim(}# #{cmd:)} {cmd:norow} {cmd:nocolumn} {cmd:xnegate} {cmd:ynegate} {cmd:origin}), {cmd:mcaplot} and {cmd:mcaprojection} [varlist] (after {cmd:mca}), {cmd:mdsconfig} and {cmd:mdsshepard} (after {cmd:mds}, {cmd:mdsmat} and {cmd:mdslong}), {cmd:procoverlay} (after {cmd:procrustes}){p_end}
{phang2}{cmd:predict} newvarlist (after {cmd:mds}: the configuration; after {cmd:procrustes}: {cmd:fitted} or {cmd:residuals}){p_end}
