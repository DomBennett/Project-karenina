# How Do Evolutionary Distinct Species Arise?

Are evolutionarily distinct species distinct for much the same reasons, or
are each evolutionarily distinct species distinct in their own unique way?
We answer this by modelling the evolutioanry distinctness of clades between
two time points.

![](https://upload.wikimedia.org/wikipedia/commons/8/8a/Porc_formiguer.JPG)

> Aardvark (*Orycteropus afer*): An evolutionarily distinct mammal featured
here to make the page more visually interesting.

## Details

**Study group**: mammals

**Data**: phylogeny and fossil

**OS**: UNIX

## Reproduce

Download repository and initial data files (`0_data`). Once directory structure is
setup, pipeline can be re-run from terminal:

```{bash}
Rscript run.R &> log &
```

## Process

* Download fossil records
* Add fossils stochastically to time-calibrated molecular phylogeny with a taxonomy constraint
* Calculate evolutionary distinctness at different different epochs for all clades
* Model change in evolutionary distinctness of clades between epochs

## Stages

* 1_pin: add fossils to molecular phylogeny
* 2_slice: calculate evolutionary distinctness at different epochs
* 3_wrngl: wrangle and merge data for modelling
* 4_model: model results

### Additional analysis

`additional_analysis` contains extra scripts for producing time slices using different simulated
birth-death trees.

## Dir Structure

```
0_data/
   mammalia.tre
stages/
   1_pin.R
   2_slice.R
   3_wrngl.R
   4_model.R
tools/
   pin_tools.R
   slice_tools.R
   wrngl_tools.R
```

Results from each stage will be saved in folders named after each stage.


## Key Packages

* [treeman](https://github.com/DomBennett/treeman)
* [paleobioDB](https://github.com/ropensci/paleobioDB)

## Data Sources

* [Bininda-Edmonds et al. (2007)](http://www.nature.com/nature/journal/v446/n7135/abs/nature05634.html)
* [The Paleobiology Database (2016)](https://paleobiodb.org/#/)

## Reference

Bennett DJ, Sutton MD, Turvey ST. 2019 How the past impacts the future: modelling the performance of evolutionarily distinct mammals through time. *Phil. Trans. R. Soc. B* 20190210. [DOI](http://dx.doi.org/10.1098/rstb.2019.0210)

## Author

D.J. Bennett

### Notes

Project name originates from the first line of [Tolstoy's "Anna Karenina"](https://en.wikipedia.org/wiki/Anna_Karenina): *Happy families are all alike; every unhappy family is unhappy in its own way.* Are evolutionarily distinct taxa each distinct in their own way?
