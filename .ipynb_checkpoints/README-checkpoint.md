# Visualizing Human Genetic Diversity

This repository contains the code used to create the blog post, [Visualizing Human Genetic Diversity](https://james-kitchens.com/blog/visualizing-human-genetic-diversity), on my personal website. This post revisualizes the dataset used in [Biddanda et al. 2020](https://doi.org/10.7554/eLife.60107), which can be found [here](https://datadryad.org/stash/dataset/doi:10.5061/dryad.rjdfn2z7v). Place the Dryad download within the assets folder and change its name to Dryad.

generate_geovar_code.py requires the `geovar` package (see this [tutorial](https://aabiddanda.github.io/geovar/index.html)) for installation instructions. This will create a new file within assets with counts for the geovar codes.

That file is an input for create_euler_diagrams.R. This code requires the `eulerr` and `UpSetR` packages. Code associated with each figure has been marked. The `eulerr` package calculates the position and shape of the ellipses within the Euler diagrams. This is not the best workflow, but these sections output a JSON object that I then copy and paste into the Markdown. Due to random starting conditions with `eulerr`, your Euler diagrams may not match those in post, but general patterns should be preserved.

visualizing-human-genetic-diversity.md is the Markdown file that is used to generate the post. The section within the "---" is the front matter, which in most scenarios can be removed or ignored. Styling of the figures may be slightly changed from the output of the R code.

If you have any questions, feel free to reach contact me (**James Kitchens**) at kitchensjn@gmail.com.