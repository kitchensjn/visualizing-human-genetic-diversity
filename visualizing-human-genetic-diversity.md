---
layout: "post"
permalink: "/blog/visualizing-human-genetic-diversity"
title: "Visualizing Human Genetic Diversity"
date: "May 16, 2023"
authors: "James Kitchens and Graham Coop"
skills: [Genetics, D3, Python, R]
description: A key insight from human genetics is that, as a species, we are all very genetically similar to one another and share much of our genetic variation. Our genome can be depicted as a string of letters (A, T, G, and C), referring to the four nucleobases found in DNA.
image: "/assets/blog/visualizing-human-genetic-diversity/thumbnail.png"
twitter:
    card: summary
---

<br>

<script src="https://d3js.org/d3.v4.js"></script>
<style>
    .euler_tooltip {
        position: absolute;
        max-width: 500px;
        padding: 10px;
        background-color: white;
        color: black;
        text-align: left;
        z-index: 9999;
    }
    path {
        transition: opacity 0.25s, stroke-width 0.25s;
    }
    svg {
        display: block;
    }
    #all_euler_diagrams {
        display: flex;
        flex-flow: row wrap;
        justify-content: space-evenly;
    }
    .upset {
        text-align: center;
    }
    .upset img {
        width: 100%;
        max-width: 1000px;
    }

    .hidden_thumbnail {
        display: none;
    }
</style>


<script>
    var stepBefore = d3.line().curve(d3.curveStepBefore);
    // Define the div for the tooltip
    var tip = d3.select(".section").append("div")
        .attr("class", "euler_tooltip")
        .style("display", "none");

    function wrap(
        text,
        width,
        dyAdjust,
        lineHeightEms,
        lineHeightSquishFactor,
        splitOnHyphen,
        centreVertically
        ) {
        // Use default values for the last three parameters if values are not provided.
        if (!lineHeightEms) lineHeightEms = 1.05;
        if (!lineHeightSquishFactor) lineHeightSquishFactor = 1;
        if (splitOnHyphen == null) splitOnHyphen = true;
        if (centreVertically == null) centreVertically = true;

        text.each(function () {
            var text = d3.select(this),
            x = text.attr("x"),
            y = text.attr("y");

            var words = [];
            text
            .text()
            .split(/\s+/)
            .forEach(function (w) {
                if (splitOnHyphen) {
                var subWords = w.split("-");
                for (var i = 0; i < subWords.length - 1; i++)
                    words.push(subWords[i] + "-");
                words.push(subWords[subWords.length - 1] + " ");
                } else {
                words.push(w + " ");
                }
            });

            text.text(null); // Empty the text element

            // `tspan` is the tspan element that is currently being added to
            var tspan = text.append("tspan");

            var line = ""; // The current value of the line
            var prevLine = ""; // The value of the line before the last word (or sub-word) was added
            var nWordsInLine = 0; // Number of words in the line
            for (var i = 0; i < words.length; i++) {
            var word = words[i];
            prevLine = line;
            line = line + word;
            ++nWordsInLine;
            tspan.text(line.trim());
            if (tspan.node().getComputedTextLength() > width && nWordsInLine > 1) {
                // The tspan is too long, and it contains more than one word.
                // Remove the last word and add it to a new tspan.
                tspan.text(prevLine.trim());
                prevLine = "";
                line = word;
                nWordsInLine = 1;
                tspan = text.append("tspan").text(word.trim());
            }
            }

            var tspans = text.selectAll("tspan");

            var h = lineHeightEms;
            // Reduce the line height a bit if there are more than 2 lines.
            if (tspans.size() > 2)
            for (var i = 0; i < tspans.size(); i++) h *= lineHeightSquishFactor;

            tspans.each(function (d, i) {
            // Calculate the y offset (dy) for each tspan so that the vertical centre
            // of the tspans roughly aligns with the text element's y position.
            var dy = i * h + dyAdjust;
            if (centreVertically) dy -= ((tspans.size() - 1) * h) / 2;
            d3.select(this)
                .attr("y", y)
                .attr("x", x)
                .attr("dy", dy + "em");
            });
        });
        }

    function d3_euler(data, id, dim, title, title_color, flip_title, loc, interactive) {

        var svg_container = d3.select(loc)
            .append("svg")
            .attr("id", id)
            .attr("width", dim)
            .attr("viewBox", "0 0 " + dim + " " + dim);

        var ellipses = data;

        var downscale = 15;
        var shift = dim/2;

        svg_container
            .selectAll("ellipses")
            .data(ellipses)
            .enter()
            .append("path")
            .attr("id", function(d) { return d._row; })
            .attr("stroke", function(d) { return d.color; })
            .attr("stroke-width", 3)
            .attr("stroke-dasharray", function(d) { return d.stroke_dasharray; })
            .attr("fill", function(d) { return d.fill; })
            .attr("d", function(d) {
                return ["M", (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift), "A", d.a/downscale, d.b/downscale, 0, 1, 1, (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift-0.001)].join(' ')
            })
            .attr("transform", function(d) {
                return "rotate(" + d.phi*(180/Math.PI) + " " + (d.h/downscale+shift) + " " + (d.k/downscale+shift) + ")"; 
            });
        
        if (interactive) {
            svg_container
                .selectAll("path")
                    .on("mouseover", function(d) {
                        d3.select("#"+id)
                            .selectAll("path")
                                .style("opacity", 0.25);
                        d3.select(this)
                            .style("cursor", "pointer")
                            .style("stroke-width", 7)
                            .style("opacity", 1);
                        tip
                            .style("display", "block")
                            .html("<h3 style='margin: 0px;'>" + d.description + "</h3></br>1KGP Abbreviation: " + d.abbreviation + "</br>Number of Sampled Individuals: " + d.sampled_individuals.toLocaleString() + "</br>Number of Common Variants: " + d.common_variants.toLocaleString() + "</br>Number of Unshared Common Variants: "+d.unshared_common_variants.toLocaleString())
                            .style("border", d.color + " solid 7px")
                            .style("left", (d3.event.pageX) + "px")
                            .style("top", (d3.event.pageY+25) + "px");
                    })
                    .on("mouseout", function(d) {
                        tip.style("display", "none");
                        d3.select("#"+id)
                            .selectAll("path")
                                .style("opacity", 1)
                                .style("stroke-width", 3);
                    });
        }
        
        if (title != "") {
            if (flip_title != "") {
                svg_container
                    .append("text")
                    .text(title)
                    .attr("x", dim/2)
                    .attr("y", dim/2)
                    .attr("text-anchor", "middle")
                    .attr("transform", "rotate(180, " + dim/2 + ", " + dim/2 + ")")
                    .style("font-family", "Arial")
                    .style("font-size", 25)
                    .style("fill", title_color)
                    .call(wrap, 300, 0, 1);
            } else {
                svg_container
                    .append("text")
                    .text(title)
                    .attr("x", dim/2)
                    .attr("y", dim/2)
                    .attr("text-anchor", "middle")
                    .style("font-family", "Arial")
                    .style("font-size", 25)
                    .style("fill", title_color)
                    .call(wrap, 300, 0, 1);
            }
        }
    }

    function add_coffee_stain_interactivity(id) {
        d3.select("#"+id)
            .selectAll("path")
                .on("mouseover", function(d) {
                    d3.select("#"+id)
                        .selectAll("path")
                            .style("cursor", "pointer")
                            .style("opacity", 0.25);
                    d3.select("#"+id)
                        .select("path")
                            .style("opacity", 1);
                })
                .on("mouseout", function(d) {
                    tip.style("display", "none");
                    d3.select("#"+id)
                        .selectAll("path")
                            .style("opacity", 1);
                });
    }    
</script>


A key insight from human genetics is that, as a species, we are all very genetically similar to one another and share much of our genetic variation. Our genome can be depicted as a string of letters (A, T, G, and C), referring to the four nucleobases found in DNA. Two human genomes picked at random are identical at ~99.9% of sites (e.g. [Mallick et al. 2016](https://www.nature.com/articles/nature18964)); in that small fraction that doesn't match (~1/1000 sites), your chromosome might carry an A while the other person’s chromosome carries a T. The majority of sites with variation have no known function; indeed, carrying an A instead of an T may have no discernible effect on your traits. Much of the common genetic variation is shared among human groups ([Lewontin, 1972](https://www.vanderbilt.edu/evolution/wp-content/uploads/sites/295/2022/04/lewontin1972.pdf)). Human geneticists are interested both in understanding which sites in the genome are functional and in unraveling the subtle differences between individuals and groups that highlight our shared history. 

Due to the complexity of data involved, understanding and visualizing patterns of human genetic variation is often challenging. One helpful place to start is to visualize the global frequencies of variants at individual sites within the genome to see how variation is shared - check out the [Geography of Genetic Variants Browser](https://popgen.uchicago.edu/ggv/) from the Novembre Lab for a nice interactive tool ([Marcus & Novembre, 2016](https://doi.org/10.1093/bioinformatics/btw643)). However, because the human genome contains approximately 3 billion sites, it would take a few lifetimes to walk through the genome  in this manner, so researchers often turn to genome-wide summary statistics to capture patterns of genetic variation.

We wanted to share some resources that we’ve been putting together for teaching human genetics using data from the 1000 Genomes Project, inspired by [Donovan et al. (2019)](https://doi.org/10.1002/sce.21506) and [Biddanda et al. (2020)](https://doi.org/10.7554/eLife.60107). These visualizations first center on the variation in a set of diverse samples from the Americas (see Figure 2) before expanding to include more globally distributed examples. In a small sample of people, we expect that they vary at only a small fraction of sites in their entire sequenced genomes.<sup>1</sup> Most of this variation is rare, and though these rare variants can be medically salient, they are the properties of specific people and their immediate families, rather than of the larger human groups. To learn about more widely shared variation and following methods similar to those in Biddanda et al., we defined a variant as “common” in a sample if it was found in more than 5% of people’s chromosomes and then filtered the data based on this criterion.

<br>
<div id="figure1" style="display: flex; justify-content: center;"></div>

<script>
    var dim = 600;

    var svg_container = d3.select("#figure1")
        .append("svg")
        .attr("id", "americas_genome")
        .attr("width", dim)
        .attr("viewBox", "0 0 " + dim + " " + dim);

    
    var ellipses = [
        {"h":-173.4314,"k":0,"a":30382.5389,"b":30382.5389,"phi":5.0914,"color":"#A9A9A9","fill":"#FFFFFF","stroke_dasharray":"5","_row":"Measurable"},
        {"h":9660.3273,"k":0,"a":3527.7784,"b":3527.7784,"phi":5.0914,"color":"#A9A9A9","fill":"#FFFFFF","stroke_dasharray":"5","_row":"Variants"},
        {"h":10754.3423,"k":362.1589,"a":1783.4049,"b":1783.4049,"phi":5.0914,"color":"none","fill":"#56B4E9","stroke_dasharray":"none","_row":"Common"}
    ] 

    var downscale = 130;
    var shift = dim/2;

    var measurable = svg_container
        .selectAll("ellipses")
        .data(ellipses)
        .enter()
        .filter(function(d) { return (d._row == "Measurable"); })
        .append("g");

    measurable
        .append("path")
        .attr("d", function(d) { return stepBefore([[d.h/downscale+shift-200, d.k/downscale+d.b/downscale+shift-480], [d.h/downscale+shift, d.k/downscale+d.b/downscale+shift-400]]) ;})
        .attr("stroke-width", 2)
        .attr("stroke", "#A9A9A9")
        .attr("fill", "none");

    measurable
        .append("text")
        .text("Measurable locations in the human genome (~2.9 billion sites)")
        .attr("x", function(d) {return d.h/downscale+shift;})
        .attr("y", function(d) {return d.k/downscale+d.b/downscale+shift-490; })
        .style("text-anchor", "middle")
        .style("fill", "#696969")
        .style("font-family", "Arial")
        .style("font-size", 18);

    measurable
        .append("path")
        .attr("id", function(d) { return d._row; })
        .attr("stroke", function(d) { return d.color; })
        .attr("stroke-width", 3)
        .attr("stroke-dasharray", function(d) { return d.stroke_dasharray; })
        .attr("fill", function(d) { return d.fill; })
        .attr("d", function(d) {
            return ["M", (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift), "A", d.a/downscale, d.b/downscale, 0, 1, 1, (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift-0.001)].join(' ')
        })
        .attr("transform", function(d) {
            return "rotate(" + d.phi*(180/Math.PI) + " " + (d.h/downscale+shift) + " " + (d.k/downscale+shift) + ")"; 
        });

   var variants = svg_container
        .selectAll("ellipses")
        .data(ellipses)
        .enter()
        .filter(function(d) { return (d._row == "Variants"); })
        .append("g");

    variants
        .append("path")
        .attr("d", function(d) { return stepBefore([[d.h/downscale+shift-100, d.k/downscale+d.b/downscale+shift-65], [d.h/downscale+shift, d.k/downscale+d.b/downscale+shift-25]]) ;})
        .attr("stroke-width", 2)
        .attr("stroke", "#A9A9A9")
        .attr("fill", "none");

    variants
        .append("text")
        .text("Observed variants within 609 individuals from the Americas (~39 million sites)")
        .attr("x", function(d) {return d.h/downscale+shift-100;})
        .attr("y", function(d) {return d.k/downscale+d.b/downscale+shift-110; })
        .style("text-anchor", "middle")
        .style("fill", "#696969")
        .style("font-family", "Arial")
        .style("font-size", 18)
        .call(wrap, 250, 1, 1);

    variants
        .append("path")
        .attr("id", function(d) { return d._row; })
        .attr("stroke", function(d) { return d.color; })
        .attr("stroke-width", 3)
        .attr("stroke-dasharray", function(d) { return d.stroke_dasharray; })
        .attr("fill", function(d) { return d.fill; })
        .attr("d", function(d) {
            return ["M", (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift), "A", d.a/downscale, d.b/downscale, 0, 1, 1, (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift-0.001)].join(' ')
        })
        .attr("transform", function(d) {
            return "rotate(" + d.phi*(180/Math.PI) + " " + (d.h/downscale+shift) + " " + (d.k/downscale+shift) + ")"; 
        });

    var common = svg_container
        .selectAll("text")
        .data(ellipses)
        .enter()
        .filter(function(d) { return (d._row == "Common"); })
        .append("g");

    common
        .append("path")
        .attr("d", function(d) { return stepBefore([[d.h/downscale+shift, d.k/downscale+d.b/downscale+shift-10], [d.h/downscale+shift-20, d.k/downscale+d.b/downscale+shift+60]]) ;})
        .attr("stroke-width", 2)
        .attr("stroke", "#56B4E9")
        .attr("fill", "none");

    common
        .append("text")
        .text("Common variants within 609 individuals from the Americas (~10 million sites)")
        .attr("x", function(d) {return d.h/downscale+shift-150;})
        .attr("y", function(d) {return d.k/downscale+d.b/downscale+shift+50; })
        .style("text-anchor", "middle")
        .style("fill", "#56B4E9")
        .style("font-family", "Arial")
        .style("font-size", 18)
        .call(wrap, 250, 1, 1);

    common
        .append("path")
        .attr("id", function(d) { return d._row; })
        .attr("stroke", function(d) { return d.color; })
        .attr("stroke-width", 3)
        .attr("stroke-dasharray", function(d) { return d.stroke_dasharray; })
        .attr("fill", function(d) { return d.fill; })
        .attr("d", function(d) {
            return ["M", (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift), "A", d.a/downscale, d.b/downscale, 0, 1, 1, (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift-0.001)].join(' ')
        })
        .attr("transform", function(d) {
            return "rotate(" + d.phi*(180/Math.PI) + " " + (d.h/downscale+shift) + " " + (d.k/downscale+shift) + ")"; 
        });

</script>

<center style="font-style: italic;">
    <p style="font-size: 16px; max-width: 1000px; text-align: left;">
        <span style="font-weight: bold;">Figure 1 - The scale of common variants in the Americas compared to the human genome.</span> The area of each circle is scaled proportionally by the number of sites in that category. The small blue circle corresponds with the number of common variants; “common” is defined as having a minor allele frequency of greater than 5% in at least one of the samples.
    </p>
</center>
<br>

The small blue circle in the above figure captures just how little variation rises to this frequency in the Americas. As the rest of this post focuses on the sharing of these common variants, it’s important to maintain perspective regarding the scale of these differences relative to the size of the human genome.

There are seven different samples from the Americas in the 1000 Genomes Project dataset (as described in Biddanda et al.), each sample being made up of 60-105 people, and we counted the number of common variants found in each sample.<sup>2</sup>

<br>
<div id="figure2" style="display: flex; justify-content: center;"></div>

<script>
    var dim = 750;

    var svg_container = d3.select("#figure2")
        .append("svg")
        .attr("id", "americas_overview")
        .attr("width", dim)
        .attr("viewBox", "0 0 " + dim + " " + dim);

    var ellipses = [{"abbreviation":"PEL","h":2493.9592,"k":-3127.3259,"a":1279.1135,"b":1279.1135,"phi":0,"common_variants":5140058,"unshared_common_variants":99338,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"#E69F00","fill":"none","stroke_dasharray":"none"},{"abbreviation":"MXL","h":2493.9592,"k":3127.3259,"a":1342.6299,"b":1342.6299,"phi":0,"common_variants":5663208,"unshared_common_variants":43322,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"},{"abbreviation":"CEU","h":-3603.8755,"k":1735.535,"a":1350.0972,"b":1350.0972,"phi":0,"common_variants":5726377,"unshared_common_variants":184313,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},{"abbreviation":"CLM","h":-3603.8755,"k":-1735.535,"a":1369.9063,"b":1369.9063,"phi":0,"common_variants":5895649,"unshared_common_variants":37361,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"#009E73","fill":"none","stroke_dasharray":"none"},{"abbreviation":"PUR","h":-890.0837,"k":3899.7116,"a":1388.1292,"b":1388.1292,"phi":0,"common_variants":6053543,"unshared_common_variants":51810,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"#000000","fill":"none","stroke_dasharray":"none"},{"abbreviation":"ASW","h":-890.0837,"k":-3899.7116,"a":1558.5375,"b":1558.5375,"phi":0,"common_variants":7631052,"unshared_common_variants":321566,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},{"abbreviation":"ACB","h":4000,"k":0,"a":1597.628,"b":1597.628,"phi":0,"common_variants":8018649,"unshared_common_variants":840969,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"#D55E00","fill":"none","stroke_dasharray":"none"}];

    /*
    var ellipses = [
        {"abbreviation":"PEL","h":2493.9592,"k":-3127.3259,"a":1279.1135,"b":1279.1135,"phi":0,"common_variants":5140058,"unshared_common_variants":99338,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"#E69F00","fill":"none","stroke_dasharray":"none"},
        {"abbreviation":"MXL","h":-890.0837,"k":-3899.7116,"a":1342.6299,"b":1342.6299,"phi":0,"common_variants":5663208,"unshared_common_variants":43322,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"},
        {"abbreviation":"CEU","h":-3603.8755,"k":-1735.535,"a":1350.0972,"b":1350.0972,"phi":0,"common_variants":5726377,"unshared_common_variants":184313,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
        {"abbreviation":"CLM","h":-3603.8755,"k":1735.535,"a":1369.9063,"b":1369.9063,"phi":0,"common_variants":5895649,"unshared_common_variants":37361,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
        {"abbreviation":"PUR","h":-890.0837,"k":3899.7116,"a":1388.1292,"b":1388.1292,"phi":0,"common_variants":6053543,"unshared_common_variants":51810,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"#000000","fill":"none","stroke_dasharray":"none"},
        {"abbreviation":"ASW","h":2493.9592,"k":3127.3259,"a":1558.5375,"b":1558.5375,"phi":0,"common_variants":7631052,"unshared_common_variants":321566,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
        {"abbreviation":"ACB","h":4000,"k":0,"a":1597.628,"b":1597.628,"phi":0,"common_variants":8018649,"unshared_common_variants":840969,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"#D55E00","fill":"none","stroke_dasharray":"none"}
    ];
    */

    var downscale = 15;
    var shift = dim/2;

    svg_container
        .append("text")
        .text("Common genetic variants found within samples from the Americas")
        .attr("x", dim/2)
        .attr("y", dim/2)
        .attr("text-anchor", "middle")
        .style("font-family", "Arial")
        .style("font-size", 25)
        .style("fill", "#696969")
        .call(wrap, 300, -1, 1)
    
    svg_container
        .append("text")
        .text('"Common" is defined as having a minor allele frequency >5%.')
        .attr("x", dim/2)
        .attr("y", dim/2)
        .attr("text-anchor", "middle")
        .style("font-family", "Arial")
        .style("font-size", 14)
        .style("fill", "#696969")
        .call(wrap, 300, 4, 1)

    var ellipse_container = svg_container
        .selectAll("ellipses")
        .data(ellipses)
        .enter()
        .append("g")
        .attr("id", function(d) { return d.abbreviation; });

    ellipse_container
        .append("path")
        .attr("stroke", function(d) { return d.color; })
        .attr("stroke-width", 3)
        .attr("stroke-dasharray", function(d) { return d.stroke_dasharray; })
        .attr("fill", function(d) { return d.fill; })
        .attr("d", function(d) {
            return ["M", (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift), "A", d.a/downscale, d.b/downscale, 0, 1, 1, (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift-0.001)].join(' ')
        })
        .attr("transform", function(d) {
            return "rotate(" + d.phi*(180/Math.PI) + " " + (d.h/downscale+shift) + " " + (d.k/downscale+shift) + ")"; 
        });

    ellipse_container
        .append("text")
        .text(function(d) { return d.description; })
        .attr("x", function(d) { return d.h/downscale+shift; })
        .attr("y", function(d) { return d.k/downscale+shift; })
        .attr("text-anchor", "middle")
        .style("fill", function(d) { return d.color; })
        .style("font-family", "Arial")
        .style("font-weight", "bold")
        .style("font-size", 14)
        .call(wrap, 150, -0.75, 1);
    
    ellipse_container
        .append("text")
        .text(function(d) { return d.sampled_individuals.toLocaleString() + " individuals"; })
        .attr("x", function(d) { return d.h/downscale+shift; })
        .attr("y", function(d) { return d.k/downscale+shift; })
        .attr("text-anchor", "middle")
        .attr("dominant-baseline", "middle")
        .style("fill", function(d) { return d.color; })
        .style("font-family", "Arial")
        .style("font-size", 14)
        .call(wrap, 150, 2, 1);
        
    ellipse_container
        .append("text")
        .text(function(d) { return d.common_variants.toLocaleString() + " variants"; })
        .attr("x", function(d) { return d.h/downscale+shift; })
        .attr("y", function(d) { return d.k/downscale+shift; })
        .attr("text-anchor", "middle")
        .attr("dominant-baseline", "middle")
        .style("fill", function(d) { return d.color; })
        .style("font-family", "Arial")
        .style("font-size", 14)
        .call(wrap, 150, 3.5, 1);
</script>

<center style="font-style: italic;">
    <p style="font-size: 16px; max-width: 1000px; text-align: left;">
        <span style="font-weight: bold;">Figure 2 - Number of common variants in seven samples from the Americas.</span> The area of each circle is proportional to the number of common variants within that sample from the 1000 Genomes Project. A “common” variant is defined as having a minor allele frequency of greater than 5%, where the minor allele identity is determined by its global allele frequency (its frequency across all samples in the 1000 Genomes Project). The number of individuals within each sample has also been included to ensure that this quantity is relatively consistent between samples.
    </p>
</center>
<br>

The levels of genetic diversity, shown as differences in the number of common variants, vary between samples: African Caribbean in Barbados (ACB) and African Ancestry in Southwest US (ASW) display the highest levels of variation. Similar to Figure 1 from Donovan et al. (2019), we implement an Euler diagram to visualize the amount of overlap in common genetic variation between samples (Figure 3). This style of visualization is like a Venn diagram, with the added property that the areas and overlaps of the shapes are proportional to the number of common variants in the corresponding samples.

<br>
<div id="figure3" style="display: flex; justify-content: center;"></div>

<script>
    d3_euler(
        data = [
            {"abbreviation":"ACB","h":-351.0253,"k":-117.9368,"a":1661.9991,"b":1482.8142,"phi":-2.3891,"common_variants":8018649,"unshared_common_variants":840969,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-186.6873,"k":-117.9368,"a":1753.6272,"b":1346.1266,"phi":-1.7434,"common_variants":7631052,"unshared_common_variants":321566,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":83.4634,"k":437.6332,"a":1518.3666,"b":1150.6006,"phi":1.187,"common_variants":5726377,"unshared_common_variants":184313,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":163.3362,"k":321.562,"a":1247.0832,"b":1433.2207,"phi":-0.2126,"common_variants":5895649,"unshared_common_variants":37361,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":188.7347,"k":343.7469,"a":1298.0515,"b":1327.9438,"phi":-1.3532,"common_variants":5663208,"unshared_common_variants":43322,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":245.9637,"k":354.8639,"a":1139.4186,"b":1374.1214,"phi":1.0056,"common_variants":5140058,"unshared_common_variants":99338,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"#E69F00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":48.8014,"k":359.7592,"a":1435.6344,"b":1296.4881,"phi":-1.5072,"common_variants":6053543,"unshared_common_variants":51810,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"#000000","fill":"none","stroke_dasharray":"none"}
        ],
        id = "americas_example_euler",
        dim = 350,
        title = "Americas",
        title_color = "#696969",
        flip_title = "",
        loc = "#figure3",
        interactive = true
    )
</script>

<center style="font-style: italic;">
    <p style="font-size: 16px; max-width: 1000px; text-align: left;">
        <span style="font-weight: bold;">Figure 3 - Overlap in common variants between samples from the Americas.</span> This Euler diagram is interactive; hovering over the ellipses displays a popup with information about the samples. It is not mathematically possible to generate ellipses with a given overlap without distortions to the areas. See the Technical details section (below) for statistics quantifying the slight errors in this and following Euler diagrams.
    </p>
</center>
<br>

It’s clear from this figure that the majority of common variants are not unique to a single sample. Instead, they are often widely distributed and shared between samples, resulting in a large degree of overlap between ellipses. The African Caribbean and African American (ACB and ASW) samples share nearly all of the common variation found in other samples. However, as noted above, they also have greater amounts of genetic variation compared to that found in the other samples (larger area), and some of that variation is not common in the other samples from the Americas. This does not mean that these variants are completely absent from the other groups, but instead, that these variants are rare or undetected in the other samples included in the figure. For example, maybe 10% of people’s chromosomes in the ACB sample carry a T instead of an A at a particular site, but this T is found in only 1% of the CEU sample.

To look at the overlap in a different way, we first considered the variation that is common (>5%) in a given sample and then identified in which other samples the variant is also common.

<br>
<div id="figure4" style="display: flex; justify-content: center;"></div>

<script>
    d3.select("#figure4")
        .append("div")
            .attr("id", "common_pops")
            .style("display", "flex")
            .style("flex-flow", "row wrap")
            .style("justify-content", "space-evenly");

    var common_dim = 300;

    d3_euler(
        data = [
            {"abbreviation":"ASW","h":-33.6008,"k":-66.9563,"a":1628.118,"b":1350.1618,"phi":-2.5773,"common_variants":7020928,"unshared_common_variants":0,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":68.6404,"k":95.4167,"a":1053.3382,"b":1381.2585,"phi":-2.9205,"common_variants":4666113,"unshared_common_variants":0,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":34.2802,"k":9.8817,"a":1448.5815,"b":1060.8529,"phi":-1.4555,"common_variants":4896250,"unshared_common_variants":0,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":16.8686,"k":59.5652,"a":1401.8114,"b":1044.7415,"phi":1.5109,"common_variants":4672565,"unshared_common_variants":0,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":-68.7641,"k":66.9836,"a":988.4922,"b":1348.208,"phi":-3.1739,"common_variants":4261270,"unshared_common_variants":0,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"#E69F00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":81.6678,"k":28.1504,"a":1428.3514,"b":1128.7544,"phi":-1.3564,"common_variants":5112528,"unshared_common_variants":0,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"#000000","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"ACB","h":-65.228,"k":-66.9563,"a":1514.9706,"b":1655.0866,"phi":-0.9991,"common_variants":8018649,"unshared_common_variants":840969,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"#D55E00","fill":"none","stroke_dasharray":"none"}
        ],
        id = "americas_acb",
        dim = common_dim,
        title = "ACB",
        title_color = "#D55E00",
        flip_title = "",
        loc = "#common_pops",
        interactive = true
    );
    
    d3_euler(
        data = [
            {"abbreviation":"ACB","h":-127.9212,"k":-16.3857,"a":1434.9282,"b":1516.3261,"phi":-3.5615,"common_variants":7020928,"unshared_common_variants":0,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":222.857,"k":-52.2642,"a":1161.4243,"b":1292.3007,"phi":-4.4373,"common_variants":4810015,"unshared_common_variants":0,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":202.7879,"k":32.4151,"a":1292.3632,"b":1225.6845,"phi":-3.1563,"common_variants":5043889,"unshared_common_variants":0,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":231.3601,"k":45.6133,"a":1272.4161,"b":1173.1753,"phi":-3.0802,"common_variants":4819661,"unshared_common_variants":0,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":271.5502,"k":108.1057,"a":1118.7258,"b":1215.4851,"phi":-1.7843,"common_variants":4398728,"unshared_common_variants":0,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"#E69F00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":144.5883,"k":9.885,"a":1369.9257,"b":1193.8532,"phi":-3.0874,"common_variants":5253401,"unshared_common_variants":0,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"#000000","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-43.4144,"k":-16.3857,"a":1572.7413,"b":1504.6041,"phi":-5.7499,"common_variants":7631052,"unshared_common_variants":321566,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"#0072B2","fill":"none","stroke_dasharray":"none"}
        ],
        id = "americas_asw",
        dim = common_dim,
        title = "ASW",
        title_color = "#0072B2",
        flip_title = "",
        loc = "#common_pops",
        interactive = true
    );

    d3_euler(
        data = [
            {"abbreviation":"ACB","h":-94.4748,"k":-70.8605,"a":1180.1411,"b":1233.1487,"phi":-0.0171,"common_variants":4666113,"unshared_common_variants":0,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-34.3686,"k":-70.8605,"a":1255.8523,"b":1193.2149,"phi":3.5755,"common_variants":4810015,"unshared_common_variants":0,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":21.3991,"k":14.8651,"a":1252.8401,"b":1313.3533,"phi":0.1056,"common_variants":5259585,"unshared_common_variants":0,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":7.914,"k":57.7595,"a":1228.54,"b":1279.6367,"phi":4.8983,"common_variants":4999087,"unshared_common_variants":0,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":58.5053,"k":91.4926,"a":1206.1001,"b":1147.0535,"phi":4.8701,"common_variants":4440037,"unshared_common_variants":0,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"#E69F00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":13.3255,"k":-21.4317,"a":1305.575,"b":1281.8573,"phi":2.1193,"common_variants":5329963,"unshared_common_variants":0,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"#000000","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":13.7257,"k":0,"a":1337.84,"b":1328.7352,"phi":4.81,"common_variants":5726377,"unshared_common_variants":184313,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"}
        ],
        id = "americas_ceu",
        dim = common_dim,
        title = "CEU",
        title_color = "#CC79A7",
        flip_title = "",
        loc = "#common_pops",
        interactive = true
    );

    d3_euler(
        data = [
            {"abbreviation":"ACB","h":-126.6347,"k":0,"a":1254.8892,"b":1196.7138,"phi":1.3928,"common_variants":4896250,"unshared_common_variants":0,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-79.4666,"k":0,"a":1266.0896,"b":1242.7369,"phi":1.5104,"common_variants":5043889,"unshared_common_variants":0,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":50.8403,"k":-84.4662,"a":1234.7334,"b":1307.0194,"phi":2.4081,"common_variants":5259585,"unshared_common_variants":0,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":53.2131,"k":49.9096,"a":1262.1653,"b":1312.3774,"phi":0.4765,"common_variants":5327378,"unshared_common_variants":0,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":78.6214,"k":145.9914,"a":1182.0006,"b":1254.7688,"phi":3.9076,"common_variants":4817623,"unshared_common_variants":0,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"#E69F00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":5.8608,"k":-25.6225,"a":1338.6864,"b":1279.3522,"phi":0.7696,"common_variants":5521792,"unshared_common_variants":0,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"#000000","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":0,"k":3.1296,"a":1322.6022,"b":1364.0593,"phi":-0.249,"common_variants":5895649,"unshared_common_variants":37361,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"#009E73","fill":"none","stroke_dasharray":"none"}
        ],
        id = "americas_clm",
        dim = common_dim,
        title = "CLM",
        title_color = "#009E73",
        flip_title = "",
        loc = "#common_pops",
        interactive = true
    );

    d3_euler(
        data = [
            {"abbreviation":"ACB","h":-169.4859,"k":-71.7374,"a":1245.6267,"b":1157.257,"phi":1.6699,"common_variants":4672565,"unshared_common_variants":0,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-121.0587,"k":-71.7374,"a":1248.688,"b":1210.0453,"phi":1.8526,"common_variants":4819661,"unshared_common_variants":0,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":0.43,"k":-136.8283,"a":1272.186,"b":1216.143,"phi":6.1252,"common_variants":4999087,"unshared_common_variants":0,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":-5.3669,"k":-60.049,"a":1259.3057,"b":1314.575,"phi":4.5561,"common_variants":5327378,"unshared_common_variants":0,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":16.8562,"k":65.9723,"a":1169.7654,"b":1296.0166,"phi":4.6401,"common_variants":4912197,"unshared_common_variants":0,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"#E69F00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":-28.2954,"k":-89.06,"a":1284.9118,"b":1266.5903,"phi":2.8881,"common_variants":5254464,"unshared_common_variants":0,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"#000000","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":0,"k":-54.4502,"a":1334.1963,"b":1304.4956,"phi":3.294,"common_variants":5663208,"unshared_common_variants":43322,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"}
        ],
        id = "americas_mxl",
        dim = common_dim,
        title = "MXL",
        title_color = "#56B4E9",
        flip_title = "",
        loc = "#common_pops",
        interactive = true
    );

    d3_euler(
        data = [
            {"abbreviation":"ACB","h":-136.4806,"k":-11.2297,"a":1165.5168,"b":1132.8495,"phi":2.3493,"common_variants":4261270,"unshared_common_variants":0,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-97.8328,"k":-11.2297,"a":1154.5489,"b":1198.3637,"phi":0.9624,"common_variants":4398728,"unshared_common_variants":0,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":0,"k":-49.1769,"a":1222.8294,"b":1134.1021,"phi":-2.6449,"common_variants":4440037,"unshared_common_variants":0,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":-22.0052,"k":28.7626,"a":1250.9645,"b":1203.7707,"phi":0.7879,"common_variants":4817623,"unshared_common_variants":0,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":-10.2769,"k":32.8326,"a":1217.7502,"b":1255.6227,"phi":0.7452,"common_variants":4912197,"unshared_common_variants":0,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":-22.9208,"k":0,"a":1170.9151,"b":1253.9321,"phi":1.873,"common_variants":4700471,"unshared_common_variants":0,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"#000000","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":-31.3731,"k":33.3752,"a":1252.0678,"b":1275.4104,"phi":-2.1762,"common_variants":5140058,"unshared_common_variants":99338,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"#E69F00","fill":"none","stroke_dasharray":"none"}
        ],
        id = "americas_pel",
        dim = common_dim,
        title = "PEL",
        title_color = "#E69F00",
        flip_title = "",
        loc = "#common_pops",
        interactive = true
    );

    d3_euler(
        data = [
            {"abbreviation":"ACB","h":-132.2327,"k":-44.0001,"a":1272.6429,"b":1225.1824,"phi":4.289,"common_variants":5112528,"unshared_common_variants":0,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-85.0325,"k":-44.0001,"a":1288.4448,"b":1267.3528,"phi":0.6853,"common_variants":5253401,"unshared_common_variants":0,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":78.2031,"k":-126.0023,"a":1231.5328,"b":1325.4885,"phi":2.384,"common_variants":5329963,"unshared_common_variants":0,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":56.6395,"k":-40.0285,"a":1310.9291,"b":1305.6904,"phi":-0.8658,"common_variants":5521792,"unshared_common_variants":0,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":87.7264,"k":0,"a":1268.0344,"b":1273.7931,"phi":4.7499,"common_variants":5254464,"unshared_common_variants":0,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":109.9055,"k":90.6142,"a":1224.4499,"b":1182.8382,"phi":-0.4921,"common_variants":4700471,"unshared_common_variants":0,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"#E69F00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":0,"k":-60.372,"a":1376.2273,"b":1338.3484,"phi":0.467,"common_variants":6053543,"unshared_common_variants":51810,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"#000000","fill":"none","stroke_dasharray":"none"}
        ],
        id = "americas_pur",
        dim = common_dim,
        title = "PUR",
        title_color = "#000000",
        flip_title = "",
        loc = "#common_pops",
        interactive = true
    );
</script>

<center style="font-style: italic;">
    <p style="font-size: 16px; max-width: 1000px; text-align: left;">
        <span style="font-weight: bold;">Figure 4 - Sharing of common variants found in each sample from the Americas.</span> Each diagram highlights a different sample, identified in the title, that was used to filter the variants down to only those that were common in the sample. The sizes of each plot are proportional to the number of variants included in the analysis (sizes are not proportional to previous figures). See Figure S1 near the bottom of this post for an alternative visualization of this figure.
    </p>
</center>
<br>

This method of filtering results in a Euler diagram where the ellipse of the highlighted sample completely encircles the other ellipses. A sample with greater numbers of common variants that are not common in other samples will show a larger disparity in size compared with the other ellipses. As before, these figures illustrate the high degree of sharing of variation among samples in the Americas. The African Caribbean in Barbados (ACB) and African Ancestry in Southwest US (ASW) samples contain the most genetic diversity, with some of this variation being shared only between those two samples. In comparison, there is somewhat less common variation (small diagram size) in the other samples and nearly all of it is shared.

Zooming back out and putting Figure 3 back onto the scale of the whole genome, the Euler diagram shrinks down to match the fraction of common variants in the genome.

<br>
<div id="figure5" style="display: flex; justify-content: center;"></div>

<script>
    var dim = 600;

    var svg_container = d3.select("#figure5")
        .append("svg")
        .attr("id", "americas_in_context")
        .attr("width", dim)
        .attr("viewBox", "0 0 " + dim + " " + dim);

    
    var ellipses = [
        {"h":-3473.3056,"k":0,"a":30379.7837,"b":30379.7837,"phi":4.6048,"color":"#A9A9A9","fill":"white","stroke_dasharray":"5","_row":"Measurable"},
        {"h":4680.7141,"k":0,"a":3503.9703,"b":3503.9703,"phi":4.6048,"color":"#A9A9A9","fill":"white","stroke_dasharray":"5","_row":"Variants"},
        {"h":5169.0583,"k":1234.5167,"a":1418.2626,"b":1728.5732,"phi":2.7036,"color":"#D55E00","fill":"none","stroke_dasharray":"none","_row":"ACB"},
        {"h":5011.8328,"k":1184.8433,"a":1803.0682,"b":1302.0182,"phi":4.7352,"color":"#0072B2","fill":"none","stroke_dasharray":"none","_row":"ASW"},
        {"h":4870.3866,"k":528.203,"a":1420.4929,"b":1219.6246,"phi":4.6236,"color":"#CC79A7","fill":"none","stroke_dasharray":"none","_row":"CEU"},
        {"h":4783.9403,"k":660.105,"a":1418.0123,"b":1267.8085,"phi":4.9378,"color":"#009E73","fill":"none","stroke_dasharray":"none","_row":"CLM"},
        {"h":4728.8264,"k":649.0093,"a":1253.3611,"b":1370.2936,"phi":7.0731,"color":"#56B4E9","fill":"none","stroke_dasharray":"none","_row":"MXL"},
        {"h":4659.5552,"k":657.2331,"a":1108.2349,"b":1415.1116,"phi":7.1649,"color":"#E69F00","fill":"none","stroke_dasharray":"none","_row":"PEL"},
        {"h":4895.4203,"k":647.1042,"a":1415.7074,"b":1301.8497,"phi":5.137,"color":"#000000","fill":"none","stroke_dasharray":"none","_row":"PUR"}
    ] 

    var downscale = 130;
    var shift = dim/2;


    var measurable = svg_container
        .selectAll("ellipses")
        .data(ellipses)
        .enter()
        .filter(function(d) { return (d._row == "Measurable"); })
        .append("g");

    measurable
        .append("path")
        .attr("d", function(d) { return stepBefore([[d.h/downscale+shift-200, d.k/downscale+d.b/downscale+shift-480], [d.h/downscale+shift, d.k/downscale+d.b/downscale+shift-400]]) ;})
        .attr("stroke-width", 2)
        .attr("stroke", "#A9A9A9")
        .attr("fill", "none");

    measurable
        .append("text")
        .text("Measurable locations in the human genome (~2.9 billion sites)")
        .attr("x", function(d) {return d.h/downscale+shift;})
        .attr("y", function(d) {return d.k/downscale+d.b/downscale+shift-490; })
        .style("text-anchor", "middle")
        .style("fill", "#696969")
        .style("font-family", "Arial")
        .style("font-size", 18);

    measurable
        .append("path")
        .attr("id", function(d) { return d._row; })
        .attr("stroke", function(d) { return d.color; })
        .attr("stroke-width", 3)
        .attr("stroke-dasharray", function(d) { return d.stroke_dasharray; })
        .attr("fill", function(d) { return d.fill; })
        .attr("d", function(d) {
            return ["M", (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift), "A", d.a/downscale, d.b/downscale, 0, 1, 1, (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift-0.001)].join(' ')
        })
        .attr("transform", function(d) {
            return "rotate(" + d.phi*(180/Math.PI) + " " + (d.h/downscale+shift) + " " + (d.k/downscale+shift) + ")"; 
        });

   var variants = svg_container
        .selectAll("ellipses")
        .data(ellipses)
        .enter()
        .filter(function(d) { return (d._row == "Variants"); })
        .append("g");

    variants
        .append("path")
        .attr("d", function(d) { return stepBefore([[d.h/downscale+shift-100, d.k/downscale+d.b/downscale+shift-65], [d.h/downscale+shift, d.k/downscale+d.b/downscale+shift-25]]) ;})
        .attr("stroke-width", 2)
        .attr("stroke", "#A9A9A9")
        .attr("fill", "none");

    variants
        .append("text")
        .text("Observed variants within 609 individuals from the Americas (~39 million sites)")
        .attr("x", function(d) {return d.h/downscale+shift-100;})
        .attr("y", function(d) {return d.k/downscale+d.b/downscale+shift-110; })
        .style("text-anchor", "middle")
        .style("fill", "#696969")
        .style("font-family", "Arial")
        .style("font-size", 18)
        .call(wrap, 250, 1, 1);

    variants
        .append("path")
        .attr("id", function(d) { return d._row; })
        .attr("stroke", function(d) { return d.color; })
        .attr("stroke-width", 3)
        .attr("stroke-dasharray", function(d) { return d.stroke_dasharray; })
        .attr("fill", function(d) { return d.fill; })
        .attr("d", function(d) {
            return ["M", (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift), "A", d.a/downscale, d.b/downscale, 0, 1, 1, (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift-0.001)].join(' ')
        })
        .attr("transform", function(d) {
            return "rotate(" + d.phi*(180/Math.PI) + " " + (d.h/downscale+shift) + " " + (d.k/downscale+shift) + ")"; 
        });

    var common_text = svg_container
        .selectAll("text")
        .data(ellipses)
        .enter()
        .filter(function(d) { return (d._row == "ACB"); })
        .append("g");
    
    common_text
        .append("path")
        .attr("d", function(d) { return stepBefore([[d.h/downscale+shift, d.k/downscale+d.b/downscale+shift], [d.h/downscale+shift-20, d.k/downscale+d.b/downscale+shift+60]]) ;})
        .attr("stroke-width", 2)
        .attr("stroke", "#A9A9A9")
        .attr("fill", "none");

    common_text
        .append("text")
        .text("Common variants within 609 individuals from the Americas (~10 million sites)")
        .attr("x", function(d) {return d.h/downscale+shift-150;})
        .attr("y", function(d) {return d.k/downscale+d.b/downscale+shift+50; })
        .style("text-anchor", "middle")
        .style("fill", "#696969")
        .style("font-family", "Arial")
        .style("font-size", 18)
        .call(wrap, 250, 1, 1);

    common_text
        .append("path")
        .attr("id", function(d) { return d._row; })
        .attr("stroke", function(d) { return d.color; })
        .attr("stroke-width", 3)
        .attr("stroke-dasharray", function(d) { return d.stroke_dasharray; })
        .attr("fill", function(d) { return d.fill; })
        .attr("d", function(d) {
            return ["M", (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift), "A", d.a/downscale, d.b/downscale, 0, 1, 1, (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift-0.001)].join(' ')
        })
        .attr("transform", function(d) {
            return "rotate(" + d.phi*(180/Math.PI) + " " + (d.h/downscale+shift) + " " + (d.k/downscale+shift) + ")"; 
        });

    var common = svg_container
        .selectAll("text")
        .data(ellipses)
        .enter()
        .filter(function(d) { return (d._row != "Measurable" & d._row != "Variants"); })
        .append("g");

    common
        .append("path")
        .attr("id", function(d) { return d._row; })
        .attr("stroke", function(d) { return d.color; })
        .attr("stroke-width", 3)
        .attr("stroke-dasharray", function(d) { return d.stroke_dasharray; })
        .attr("fill", function(d) { return d.fill; })
        .attr("d", function(d) {
            return ["M", (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift), "A", d.a/downscale, d.b/downscale, 0, 1, 1, (d.h/downscale+d.a/downscale+shift), (d.k/downscale+shift-0.001)].join(' ')
        })
        .attr("transform", function(d) {
            return "rotate(" + d.phi*(180/Math.PI) + " " + (d.h/downscale+shift) + " " + (d.k/downscale+shift) + ")"; 
        });
</script>

<center style="font-style: italic;">
    <p style="font-size: 16px; max-width: 1000px; text-align: left;">
        <span style="font-weight: bold;">Figure 5 - Common variants in perspective.</span> An Euler diagram of the common variants in samples located in the Americas relative to the scale of the human genome. As a small note, the positions and orientations of ellipses within the Euler diagram differ slightly from Figure 3. This is because the eulerr package gives varied results with each run due to random starting conditions within the algorithm.
    </p>
</center>
<br>

Genetic diversity in the Americas reflects the history of colonialism and the transatlantic slave trade, which has moved people from across the globe into the region over the past few hundred years. Given this, you may wonder whether the high degree of overlap reflects this recent history of the Americas or whether it is representative of sharing that is present in geographically distant samples. To look into this question, we created a Euler diagram with five samples, one from each of the broad geographic groupings used by Biddanda et al.

<br>
<div id="figure6" style="display: flex; justify-content: center; transform: rotate(180deg);"></div>

<script>
    d3_euler(
        data = [
            {"abbreviation":"BEB","h":-146.4612,"k":-157.1359,"a":1275.6961,"b":1409.5007,"phi":6.0493,"common_variants":5817047,"unshared_common_variants":184032,"description":"Bengali in Bangladesh","sampled_individuals":86,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CHB","h":-264.8249,"k":-157.1359,"a":1458.5458,"b":1072.357,"phi":5.822,"common_variants":5214970,"unshared_common_variants":251564,"description":"Han Chinese in Beijing, China","sampled_individuals":103,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"GBR","h":102.6356,"k":-237.4416,"a":1387.917,"b":1243.9952,"phi":8.025,"common_variants":5592182,"unshared_common_variants":249975,"description":"British in England and Scotland","sampled_individuals":91,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":-15.8317,"k":-188.1657,"a":1323.2218,"b":1260.246,"phi":8.806,"common_variants":5663208,"unshared_common_variants":167945,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"YRI","h":208.4168,"k":578.8131,"a":1470.648,"b":1705.8368,"phi":8.7792,"common_variants":8138465,"unshared_common_variants":3120507,"description":"Yoruba in Ibadan, Nigeria","sampled_individuals":108,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"}
        ],
        id = "global_samples",
        dim = 350,
        title = "Global",
        title_color = "#696969",
        flip_title = true,
        loc = "#figure6",
        interactive = true
    )
</script>

<center style="font-style: italic;">
    <p style="font-size: 16px; max-width: 1000px; text-align: left;">
        <span style="font-weight: bold;">Figure 6 - Overlap in common variants between samples from global samples.</span> An interactive Euler diagram of the common variants in five geographically distant samples: Bengali in Bangladesh (BEB), Han Chinese in Beijing, China (CHB), British in England and Scotland (GBR), Mexican Ancestry in Los Angeles, California (MXL), and Yoruba in Ibadan, Nigeria (YRI).
    </p>
</center>
<br>

Overall, this diagram has a very similar structure to the diagram created with the samples from the Americas. There is a high degree of overlap between all of the samples, with the higher genetic diversity of the Yoruba in Ibadan, Nigeria sample resulting in a larger ellipse that stretches outside of the cluster of other ellipses. This pattern matches the one of high diversity in the African Caribbean and African American (ACB and ASW) samples from the Americas described above. Even when considering quite geographically distant samples of humans, the dominant pattern is that of shared genetic variation. 

Lastly, given this global view, we can zoom in and look at how variation is partitioned at finer geographic scales by using all 26 samples within the 1000 Genomes Project dataset. We see that samples from Africa contain the greatest amount of genetic diversity. Much of that common genetic variation is shared, but each sample contains some variation not found in other samples. There’s a slight reduction in the variation present in samples whose recent ancestors lived outside Africa, consistent with the view that humans evolved in Africa, and when humans first migrated out of Africa, they took with them only a subset of the genetic diversity present in Africa. 

<br>
<div id="figure7" style="display: flex; justify-content: center;"></div>

<script>
    var global_dim = 300;

    d3.select("#figure7")
        .append("div")
            .attr("id", "global")
            .style("display", "flex")
            .style("flex-flow", "row wrap")
            .style("justify-content", "space-evenly");

    d3_euler(
        data = [
            {"abbreviation":"ESN","h":0,"k":0,"a":1676.94,"b":1524.893,"phi":0.0129,"common_variants":8203941,"unshared_common_variants":246520,"description":"Esan in Nigeria","sampled_individuals":99,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"GWD","h":258.3767,"k":0,"a":1451.0148,"b":1652.9013,"phi":1.4865,"common_variants":7943922,"unshared_common_variants":250670,"description":"Gambian in Western Division, The Gambia - Mandinka","sampled_individuals":113,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"LWK","h":70.8708,"k":236.3707,"a":1521.5312,"b":1659.0083,"phi":2.3397,"common_variants":8130247,"unshared_common_variants":477047,"description":"Luhya in Webuye, Kenya","sampled_individuals":99,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"MSL","h":170.9409,"k":-48.9127,"a":1772.8305,"b":1415.2643,"phi":2.1668,"common_variants":8158792,"unshared_common_variants":302274,"description":"Mende in Sierra Leone","sampled_individuals":85,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"YRI","h":40.5704,"k":-34.279,"a":1497.0224,"b":1673.7895,"phi":2.8898,"common_variants":8138465,"unshared_common_variants":168925,"description":"Yoruba in Ibadan, Nigeria","sampled_individuals":108,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"}
        ],
        id = "africa",
        dim = global_dim,
        title = "Africa",
        title_color = "#696969",
        flip_title = "",
        loc = "#global",
        interactive = true
    );

    d3_euler(
        data = [
            {"abbreviation":"FIN","h":-155.7146,"k":-108.2774,"a":1337.9647,"b":1363.3799,"phi":3.6501,"common_variants":5732641,"unshared_common_variants":281951,"description":"Finnish in Finland","sampled_individuals":99,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"GBR","h":3.7505,"k":-108.2774,"a":1316.9096,"b":1351.2426,"phi":3.2532,"common_variants":5592182,"unshared_common_variants":78046,"description":"British in England and Scotland","sampled_individuals":91,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"IBS","h":-22.1454,"k":0,"a":1389.7241,"b":1311.9542,"phi":0.737,"common_variants":5729736,"unshared_common_variants":113070,"description":"Iberian populations in Spain","sampled_individuals":107,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"TSI","h":-34.1761,"k":16.549,"a":1275.2029,"b":1429.5177,"phi":6.6987,"common_variants":5728264,"unshared_common_variants":135201,"description":"Toscani in Italy","sampled_individuals":107,"color":"#009E73","fill":"none","stroke_dasharray":"none"}
        ],
        id = "europe",
        dim = global_dim,
        title = "Europe",
        title_color = "#696969",
        flip_title = "",
        loc = "#global",
        interactive = true
    );

    d3_euler(
        data = [
            {"abbreviation":"BEB","h":0,"k":24.8724,"a":1378.1414,"b":1295.6314,"phi":3.0176,"common_variants":5817047,"unshared_common_variants":102477,"description":"Bengali in Bangladesh","sampled_individuals":86,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"GIH","h":127.3121,"k":24.8724,"a":1310.2736,"b":1363.1123,"phi":1.7942,"common_variants":5788635,"unshared_common_variants":114785,"description":"Gujarati Indians in Houston, TX","sampled_individuals":103,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"ITU","h":75.4515,"k":80.9334,"a":1360.1288,"b":1290.1062,"phi":-0.5511,"common_variants":5768108,"unshared_common_variants":72231,"description":"Indian Telugu in the UK","sampled_individuals":102,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PJL","h":57.7924,"k":0,"a":1374.8951,"b":1328.0929,"phi":1.65,"common_variants":5850455,"unshared_common_variants":120243,"description":"Punjabi in Lahore, Pakistan","sampled_individuals":96,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"STU","h":35.3873,"k":102.35,"a":1326.2806,"b":1348.3686,"phi":-0.9169,"common_variants":5771145,"unshared_common_variants":73815,"description":"Sri Lankan Tamil in the UK","sampled_individuals":102,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"}
        ],
        id = "south_asia",
        dim = global_dim,
        title = "South Asia",
        title_color = "#696969",
        flip_title = "",
        loc = "#global",
        interactive = true
    );

    d3_euler(
        data = [
            {"abbreviation":"CDX","h":-106.5464,"k":9.246,"a":1251.3989,"b":1252.4096,"phi":8.0461,"common_variants":5186192,"unshared_common_variants":80552,"description":"Chinese Dai in Xishuangbanna, China","sampled_individuals":93,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CHB","h":0,"k":9.246,"a":1292.5497,"b":1239.6108,"phi":3.0526,"common_variants":5214970,"unshared_common_variants":63247,"description":"Han Chinese in Beijing, China","sampled_individuals":103,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CHS","h":-39.4333,"k":0,"a":1304.299,"b":1253.6959,"phi":7.9788,"common_variants":5220632,"unshared_common_variants":44981,"description":"Han Chinese South","sampled_individuals":105,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"JPT","h":-19.9002,"k":117.948,"a":1293.7822,"b":1251.9766,"phi":8.3721,"common_variants":5223809,"unshared_common_variants":180251,"description":"Japanese in Tokyo, Japan","sampled_individuals":104,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"KHV","h":-107.7766,"k":6.7343,"a":1271.2812,"b":1301.8932,"phi":4.46,"common_variants":5299958,"unshared_common_variants":91705,"description":"Kinh in Ho Chi Minh City, Vietnam","sampled_individuals":99,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"}
        ],
        id = "east_asia",
        dim = global_dim,
        title = "East Asia",
        title_color = "#696969",
        flip_title = "",
        loc = "#global",
        interactive = true
    );

    d3_euler(
        data = [
            {"abbreviation":"ACB","h":-351.0253,"k":-117.9368,"a":1661.9991,"b":1482.8142,"phi":-2.3891,"common_variants":8018649,"unshared_common_variants":840969,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"#D55E00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-186.6873,"k":-117.9368,"a":1753.6272,"b":1346.1266,"phi":-1.7434,"common_variants":7631052,"unshared_common_variants":321566,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"#0072B2","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":83.4634,"k":437.6332,"a":1518.3666,"b":1150.6006,"phi":1.187,"common_variants":5726377,"unshared_common_variants":184313,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"#CC79A7","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":163.3362,"k":321.562,"a":1247.0832,"b":1433.2207,"phi":-0.2126,"common_variants":5895649,"unshared_common_variants":37361,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"#009E73","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":188.7347,"k":343.7469,"a":1298.0515,"b":1327.9438,"phi":-1.3532,"common_variants":5663208,"unshared_common_variants":43322,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"#56B4E9","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":245.9637,"k":354.8639,"a":1139.4186,"b":1374.1214,"phi":1.0056,"common_variants":5140058,"unshared_common_variants":99338,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"#E69F00","fill":"none","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":48.8014,"k":359.7592,"a":1435.6344,"b":1296.4881,"phi":-1.5072,"common_variants":6053543,"unshared_common_variants":51810,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"#000000","fill":"none","stroke_dasharray":"none"}
        ],
        id = "americas",
        dim = global_dim,
        title = "Americas",
        title_color = "#696969",
        flip_title = "",
        loc = "#global",
        interactive = true
    );

</script>

<center style="font-style: italic;">
    <p style="font-size: 16px; max-width: 1000px; text-align: left;">
        <span style="font-weight: bold;">Figure 7 - Sharing of common variation within geographic regions.</span> Five interactive Euler diagrams of the 26 global samples using the broad geographic groupings from Biddanda et al. 2020.
    </p>
</center>
<br>

It's easy for us to fall into the trap of thinking that humans are very genetically different. Historically, our ideas about the structure of human biological variation have been shaped by a few visible physical traits, notably skin color, that have a geographic pattern. But the genetic variants contributing to skin pigmentation are unrepresentative of the more general patterns of genetic sharing present among groups of people sampled from across the world. The genetic changes involved in skin pigmentation differences can show striking geographic patterns (e.g. [SLC24A5](https://popgen.uchicago.edu/ggv/?data=%221000genomes%22&chr=15&pos=48426484)), but that is because they have been shaped by strong local adaptation to the climatic conditions that people encountered as they moved around the world. These loci are fascinating examples of adaptation but are also the exception in comparison to the high degree of sharing that we see for most of human genetic variation. 

### Acknowledgements

We thank John Novembre, Arjun Biddanda, Molly Przeworski, Doc Edge, and the Coop Lab for their help during the editing of this post.


### Technical details

We used the `geovar` package in Python to group the ~92 million variants included in the 1000 Genomes Project based on minor allele frequency (MAF). Variants were separated into five bins based (MAF=0%, 0%<MAF<1%, 1%<MAF<5%, 5%<MAF<10%, and MAF>10%), though two bins would have sufficed for this analysis (MAF<5% and MAF>5%). We used the `eulerr` package in R to calculate the position and orientation of ellipses in the Euler diagrams. Unfortunately, exactly proportionally scaling the area of every region of this diagram becomes difficult to impossible as you increase the number of sets, or samples. Because of this, we have included two goodness-of-fit measurements provided by the `eulerr` package and described in further details in the package's [tutorial](https://cran.r-project.org/web/packages/eulerr/vignettes/introduction.html). For both measurements, values closer to zero have less error. Below is a table with these measurements for all of the Euler diagrams presented in this post:

<br>

<style>
    table {
        margin-left: auto;
        margin-right: auto;
        width: 100%;
        max-width: 600px;
    }

    th {
        background-color: lightgray;
        padding: 10px;
        border: solid 1px gray;
    }

    td {
        padding: 5px 10px;
        border: solid 1px gray;
    }
</style>



| figure       | stress   | diagError    |
| ------------ | -------- | ------------ |
| 3            | 4.53E-04 | 0.019573681  |
| 4_ACB        | 1.09E-04 | 0.008677686  |
| 4_ASW        | 1.52E-04 | 0.0133913503 |
| 4_CEU        | 1.04E-04 | 0.016955049  |
| 4_CLM        | 3.08E-04 | 0.0260526896 |
| 4_MXL        | 2.53E-04 | 0.0236936357 |
| 4_PEL        | 1.89E-04 | 0.0181390286 |
| 4_PUR        | 3.58E-04 | 0.0288044056 |
| 5            | 1.09E-09 | 0.0001789465 |
| 6            | 2.30E-03 | 0.0198055721 |
| 7_Africa     | 1.19E-03 | 0.027619266  |
| 7_Europe     | 6.39E-07 | 0.0004395497 |
| 7_South_Asia | 6.41E-04 | 0.0458260097 |
| 7_East_Asia  | 5.77E-04 | 0.0410841878 |
| 7_Americas   | 4.53E-04 | 0.019573681  |


<br>

The package also breaks down error by set overlap to better understand exactly which sections are over-/underrepresented by the visualization, though that is not included here. With all of that being said, these diagrams offer a unique visualization method that can be particularly useful for more qualitative interpretations of the population relationships. We converted the output of `eulerr` into a JSON format and passed this to JavaScript for plotting using D3.js. Plotting is possible directly from R, but we used D3.js for its customizability and support of interactive figures. All of the figures (alongside the code we used to generate them) can be found [here](https://github.com/kitchensjn/visualizing-human-genetic-diversity).


### Additional figures

The following figures offer alternative methods of visualization to those within this post. Details about these figures are provided in the figure captions.

<br>
<div id="figureS1" style="display: flex; justify-content: center;"></div>

<script>
    d3.select("#figureS1")
        .append("div")
            .attr("id", "coffee_stains")
            .style("display", "flex")
            .style("flex-flow", "row wrap")
            .style("justify-content", "space-evenly");

    d3_euler(
        data = [
            {"abbreviation":"ACB","h":-65.228,"k":-66.9563,"a":1514.9706,"b":1655.0866,"phi":-0.9991,"common_variants":8018649,"unshared_common_variants":840969,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"none","fill":"#D55E00","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-33.6008,"k":-66.9563,"a":1628.118,"b":1350.1618,"phi":-2.5773,"common_variants":7020928,"unshared_common_variants":0,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":68.6404,"k":95.4167,"a":1053.3382,"b":1381.2585,"phi":-2.9205,"common_variants":4666113,"unshared_common_variants":0,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":34.2802,"k":9.8817,"a":1448.5815,"b":1060.8529,"phi":-1.4555,"common_variants":4896250,"unshared_common_variants":0,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":16.8686,"k":59.5652,"a":1401.8114,"b":1044.7415,"phi":1.5109,"common_variants":4672565,"unshared_common_variants":0,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":-68.7641,"k":66.9836,"a":988.4922,"b":1348.208,"phi":-3.1739,"common_variants":4261270,"unshared_common_variants":0,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":81.6678,"k":28.1504,"a":1428.3514,"b":1128.7544,"phi":-1.3564,"common_variants":5112528,"unshared_common_variants":0,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"}
        ],
        id = "americas_acb_coffee",
        dim = common_dim,
        title = "ACB",
        title_color = "#D55E00",
        flip_title = "",
        loc = "#coffee_stains",
        interactive = false
    );
    add_coffee_stain_interactivity("americas_acb_coffee");

    d3_euler(
        data = [
            {"abbreviation":"ASW","h":-43.4144,"k":-16.3857,"a":1572.7413,"b":1504.6041,"phi":-5.7499,"common_variants":7631052,"unshared_common_variants":321566,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"none","fill":"#0072B2","stroke_dasharray":"none"},
            {"abbreviation":"ACB","h":-127.9212,"k":-16.3857,"a":1434.9282,"b":1516.3261,"phi":-3.5615,"common_variants":7020928,"unshared_common_variants":0,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":222.857,"k":-52.2642,"a":1161.4243,"b":1292.3007,"phi":-4.4373,"common_variants":4810015,"unshared_common_variants":0,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":202.7879,"k":32.4151,"a":1292.3632,"b":1225.6845,"phi":-3.1563,"common_variants":5043889,"unshared_common_variants":0,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":231.3601,"k":45.6133,"a":1272.4161,"b":1173.1753,"phi":-3.0802,"common_variants":4819661,"unshared_common_variants":0,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":271.5502,"k":108.1057,"a":1118.7258,"b":1215.4851,"phi":-1.7843,"common_variants":4398728,"unshared_common_variants":0,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":144.5883,"k":9.885,"a":1369.9257,"b":1193.8532,"phi":-3.0874,"common_variants":5253401,"unshared_common_variants":0,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"}
        ],
        id = "americas_asw_coffee",
        dim = common_dim,
        title = "ASW",
        title_color = "#0072B2",
        flip_title = "",
        loc = "#coffee_stains",
        interactive = false
    );
    add_coffee_stain_interactivity("americas_asw_coffee");

    d3_euler(
        data = [
            {"abbreviation":"CEU","h":-11.6367,"k":77.6556,"a":1308.9382,"b":1358.2787,"phi":2.7786,"common_variants":5726377,"unshared_common_variants":184313,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"none","fill":"#CC79A7","stroke_dasharray":"none"},
            {"abbreviation":"ACB","h":-117.9052,"k":0,"a":1246.4673,"b":1167.6055,"phi":7.6663,"common_variants":4666113,"unshared_common_variants":0,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-58.0648,"k":0,"a":1263.5511,"b":1186.0482,"phi":3.7958,"common_variants":4810015,"unshared_common_variants":0,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":0,"k":92.4295,"a":1233.0207,"b":1334.5812,"phi":6.0989,"common_variants":5259585,"unshared_common_variants":0,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":-25.588,"k":136.0689,"a":1266.7121,"b":1241.1807,"phi":4.0215,"common_variants":4999087,"unshared_common_variants":0,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":42.8027,"k":156.219,"a":1110.1403,"b":1246.2762,"phi":6.3725,"common_variants":4440037,"unshared_common_variants":0,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":-13.8694,"k":58.9816,"a":1270.9347,"b":1316.9558,"phi":6.1031,"common_variants":5329963,"unshared_common_variants":0,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"}
        ],
        id = "americas_ceu_coffee",
        dim = common_dim,
        title = "CEU",
        title_color = "#CC79A7",
        flip_title = "",
        loc = "#coffee_stains",
        interactive = false
    );
    add_coffee_stain_interactivity("americas_ceu_coffee");

    d3_euler(
        data = [
            {"abbreviation":"CLM","h":0,"k":3.1296,"a":1322.6022,"b":1364.0593,"phi":-0.249,"common_variants":5895649,"unshared_common_variants":37361,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"none","fill":"#009E73","stroke_dasharray":"none"},
            {"abbreviation":"ACB","h":-126.6347,"k":0,"a":1254.8892,"b":1196.7138,"phi":1.3928,"common_variants":4896250,"unshared_common_variants":0,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-79.4666,"k":0,"a":1266.0896,"b":1242.7369,"phi":1.5104,"common_variants":5043889,"unshared_common_variants":0,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":50.8403,"k":-84.4662,"a":1234.7334,"b":1307.0194,"phi":2.4081,"common_variants":5259585,"unshared_common_variants":0,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":53.2131,"k":49.9096,"a":1262.1653,"b":1312.3774,"phi":0.4765,"common_variants":5327378,"unshared_common_variants":0,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":78.6214,"k":145.9914,"a":1182.0006,"b":1254.7688,"phi":3.9076,"common_variants":4817623,"unshared_common_variants":0,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":5.8608,"k":-25.6225,"a":1338.6864,"b":1279.3522,"phi":0.7696,"common_variants":5521792,"unshared_common_variants":0,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"}
        ],
        id = "americas_clm_coffee",
        dim = common_dim,
        title = "CLM",
        title_color = "#009E73",
        flip_title = "",
        loc = "#coffee_stains",
        interactive = false
    );
    add_coffee_stain_interactivity("americas_clm_coffee");

    d3_euler(
        data = [
            {"abbreviation":"MXL","h":0,"k":-54.4502,"a":1334.1963,"b":1304.4956,"phi":3.294,"common_variants":5663208,"unshared_common_variants":43322,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"none","fill":"#56B4E9","stroke_dasharray":"none"},
            {"abbreviation":"ACB","h":-169.4859,"k":-71.7374,"a":1245.6267,"b":1157.257,"phi":1.6699,"common_variants":4672565,"unshared_common_variants":0,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-121.0587,"k":-71.7374,"a":1248.688,"b":1210.0453,"phi":1.8526,"common_variants":4819661,"unshared_common_variants":0,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":0.43,"k":-136.8283,"a":1272.186,"b":1216.143,"phi":6.1252,"common_variants":4999087,"unshared_common_variants":0,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":-5.3669,"k":-60.049,"a":1259.3057,"b":1314.575,"phi":4.5561,"common_variants":5327378,"unshared_common_variants":0,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":16.8562,"k":65.9723,"a":1169.7654,"b":1296.0166,"phi":4.6401,"common_variants":4912197,"unshared_common_variants":0,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":-28.2954,"k":-89.06,"a":1284.9118,"b":1266.5903,"phi":2.8881,"common_variants":5254464,"unshared_common_variants":0,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"}
        ],
        id = "americas_mxl_coffee",
        dim = common_dim,
        title = "MXL",
        title_color = "#56B4E9",
        flip_title = "",
        loc = "#coffee_stains",
        interactive = false
    );
    add_coffee_stain_interactivity("americas_mxl_coffee");

    d3_euler(
        data = [
            {"abbreviation":"PEL","h":-31.3731,"k":33.3752,"a":1252.0678,"b":1275.4104,"phi":-2.1762,"common_variants":5140058,"unshared_common_variants":99338,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"none","fill":"#E69F00","stroke_dasharray":"none"},
            {"abbreviation":"ACB","h":-136.4806,"k":-11.2297,"a":1165.5168,"b":1132.8495,"phi":2.3493,"common_variants":4261270,"unshared_common_variants":0,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-97.8328,"k":-11.2297,"a":1154.5489,"b":1198.3637,"phi":0.9624,"common_variants":4398728,"unshared_common_variants":0,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":0,"k":-49.1769,"a":1222.8294,"b":1134.1021,"phi":-2.6449,"common_variants":4440037,"unshared_common_variants":0,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":-22.0052,"k":28.7626,"a":1250.9645,"b":1203.7707,"phi":0.7879,"common_variants":4817623,"unshared_common_variants":0,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":-10.2769,"k":32.8326,"a":1217.7502,"b":1255.6227,"phi":0.7452,"common_variants":4912197,"unshared_common_variants":0,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"PUR","h":-22.9208,"k":0,"a":1170.9151,"b":1253.9321,"phi":1.873,"common_variants":4700471,"unshared_common_variants":0,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"}
        ],
        id = "americas_pel_coffee",
        dim = common_dim,
        title = "PEL",
        title_color = "#E69F00",
        flip_title = "",
        loc = "#coffee_stains",
        interactive = false
    );
    add_coffee_stain_interactivity("americas_pel_coffee");

    d3_euler(
        data = [
        {"abbreviation":"PUR","h":0,"k":-60.372,"a":1376.2273,"b":1338.3484,"phi":0.467,"common_variants":6053543,"unshared_common_variants":51810,"description":"Puerto Rican in Puerto Rico","sampled_individuals":104,"color":"none","fill":"#000000","stroke_dasharray":"none"},
            {"abbreviation":"ACB","h":-132.2327,"k":-44.0001,"a":1272.6429,"b":1225.1824,"phi":4.289,"common_variants":5112528,"unshared_common_variants":0,"description":"African Caribbean in Barbados","sampled_individuals":96,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"ASW","h":-85.0325,"k":-44.0001,"a":1288.4448,"b":1267.3528,"phi":0.6853,"common_variants":5253401,"unshared_common_variants":0,"description":"African Ancestry in Southwest US","sampled_individuals":61,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"CEU","h":78.2031,"k":-126.0023,"a":1231.5328,"b":1325.4885,"phi":2.384,"common_variants":5329963,"unshared_common_variants":0,"description":"Utah residents (CEPH) with Northern and Western European ancestry","sampled_individuals":99,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"CLM","h":56.6395,"k":-40.0285,"a":1310.9291,"b":1305.6904,"phi":-0.8658,"common_variants":5521792,"unshared_common_variants":0,"description":"Colombian in Medellin, Colombia","sampled_individuals":94,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"MXL","h":87.7264,"k":0,"a":1268.0344,"b":1273.7931,"phi":4.7499,"common_variants":5254464,"unshared_common_variants":0,"description":"Mexican Ancestry in Los Angeles, California","sampled_individuals":64,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"},
            {"abbreviation":"PEL","h":109.9055,"k":90.6142,"a":1224.4499,"b":1182.8382,"phi":-0.4921,"common_variants":4700471,"unshared_common_variants":0,"description":"Peruvian in Lima, Peru","sampled_individuals":85,"color":"none","fill":"#FFFFFF","stroke_dasharray":"none"}
        ],
        id = "americas_pur_coffee",
        dim = common_dim,
        title = "PUR",
        title_color = "#000000",
        flip_title = "",
        loc = "#coffee_stains",
        interactive = false
    );
    add_coffee_stain_interactivity("americas_pur_coffee");
</script>

<center style="font-style: italic;">
    <p style="font-size: 16px; max-width: 1000px; text-align: left;">
        <span style="font-weight: bold;">Figure S1 - Sharing of common variants found in each sample from the Americas.</span> Seven interactive "coffee stain" diagrams, an alternate visualization of Figure 4. The colored area is proportional in size to the number of common variants within the highlighted sample, identified in the title, that aren’t shared with another sample. Within each subfigure, the ellipse on the bottom corresponds with the highlighted sample and is filled in with that sample’s respective color. All other ellipses are filled in with white and stacked on top, thus giving the appearance of cutting out the area and leaving only the common variants that aren’t shared with another sample. If you hover over the diagram, the exact orientation of the other ellipses become more apparent.
    </p>
</center>
<br>

<br>

{:.upset}
![UpSet Plot of Five Global Samples](/assets/blog/visualizing-human-genetic-diversity/upset_global.png)

<center style="font-style: italic;">
    <p style="font-size: 16px; max-width: 1000px; text-align: left;">
        <span style="font-weight: bold;">Figure S2 - Sharing of common variation within geographic regions.</span> An UpSet plot, an alternative visualization of Figure 6. UpSet plots, created by <a href="https://doi.org/10.1109/TVCG.2014.2346248">Lex et al. (2014)</a>, are useful for handling large numbers of sets. They can communicate the exact overlap between sets, unlike Euler diagrams (as discussed in the Technical details section), but are also a bit more challenging to read as there are multiple subfigures. To draw comparisons with the Euler diagrams, the horizontal bar graph on the bottom left depicts the areas of the ellipses and the vertical bar graph shows the areas of the overlaps between ellipses referenced usings dots in the bottom subfigure.
    </p>
</center>
<br>

### Citations

Mallick et al. (2016) [The Simons Genome Diversity Project: 300 genomes from 142 diverse populations](https://doi.org/10.1038/nature18964). *Nature* 538, 201–206

Lewontin (1972) [The Apportionment of Human Diversity](https://www.vanderbilt.edu/evolution/wp-content/uploads/sites/295/2022/04/lewontin1972.pdf). Dobzhansky, T., Hecht, M.K., Steere, W.C. (eds), *Evolutionary Biology*

Marcus & Novembre (2016) [Visualizing the Geography of Genetic Variants](https://doi.org/10.1093/bioinformatics/btw643). *Bioinformatics*, Volume 33, Issue 4

Donovan et al. (2019) [Toward a more humane genetics education: Learning about the social and quantitative complexities of human genetic variation research could reduce racial bias in adolescent and adult populations](https://doi.org/10.1002/sce.21506). *Science Education*, Volume 103, Issue 3

Arjun Biddanda, Daniel P Rice, John Novembre (2020) [A variant-centric perspective on geographic patterns of human allele frequency variation](https://doi.org/10.7554/eLife.60107). *eLife*, 9:e60107

Lex et al. (2014) [UpSet: Visualization of Intersecting Sets](https://doi.org/10.1109/TVCG.2014.2346248). *IEEE Transactions on Visualization and Computer Graphics*, Volume 20, Issue 12

---

<small>
    <sup>1</sup>Indeed, if we sequence the entire population of the world, we'd see nearly every site being variable in some one. But these variants would be vanishingly rare in the population, overall.
</small>

<small>
    <sup>2</sup>It would be interesting to explore rarefaction approaches to account for the differences in the sample size ([Cotter et al. 2023](https://doi.org/10.1093/genetics/iyad070)).
</small>