# Visualizing Human Genetic Diversity

This repository contains the code used to create the blog post, [Visualizing Human Genetic Diversity](https://james-kitchens.com/blog/visualizing-human-genetic-diversity), on my personal website. This post revisualizes the dataset used in [Biddanda et al. 2020](https://doi.org/10.7554/eLife.60107), which can be found [here](https://datadryad.org/stash/dataset/doi:10.5061/dryad.rjdfn2z7v). Place the Dryad download within the assets folder and change its name to Dryad.

generate_geovar_code.py requires the `geovar` package (see this [tutorial](https://aabiddanda.github.io/geovar/index.html)) for installation instructions. This will create a new file within assets with counts for the geovar codes.

That file is an input for create_euler_diagrams.R. This code requires the `eulerr` and `UpSetR` packages. Code associated with each figure has been marked. The `eulerr` package calculates the position and shape of the ellipses within the Euler diagrams. This is not the best workflow, but these sections output a JSON object that I then copy and paste into the Markdown. Due to random starting conditions with `eulerr`, your Euler diagrams may not match those in post, but general patterns should be preserved.

visualizing-human-genetic-diversity.md is the Markdown file that is used to generate the post. The section within the "---" is the front matter, which in most scenarios can be removed or ignored. Styling of the figures may be slightly changed from the output of the R code.

## Online Activity | [Demo](https://james-kitchens.com/visualizing-human-genetic-diversity/)

We've developed a short online assignment to introduce students to Euler diagrams and their application to genetic data. If you would like to record your students' responses as part of a larger class discussion, follow the below instructions:

1) If you haven't yet, fork this repository on GitHub
2) Set up a Google Sheet and Google Apps Script to receive your students' responses:
    - [General instructions](https://github.com/jamiewilson/form-to-google-sheets) for set up
    - Add these column names to Row 1: `timestamp`, `time_spent`, `age`, `student`, `finished_codeword`, `interpretability_overlap0`, `interpretability_overlap1`, `build_your_own_overlap0`, `build_your_own_overlap1`, `thousand_genomes_overlap0`, `thousand_genomes_overlap1`, `thousand_genomes_overlap2`
    - Paste the link to your Google Sheet in Line 584 of `docs/_layouts/home.html`
    - Commit these changes to your GitHub fork
3) This activity is created with Jekyll and can be hosted through GitHub Pages:
    - Within the GitHub repository, go to `Settings` and then `Pages` in the left menu bar
    - In Build and deployment, make sure that `Source` is set to `Deploy from a branch` and `Branch` is set to `docs/`
    - In a minute, the activity should deploy. The link will be at the top of your current page.
    - Make sure to test that the answers are being recorded to your Google Sheet!

If you have any questions, feel free to reach contact me (**James Kitchens**) at kitchensjn at gmail.com.