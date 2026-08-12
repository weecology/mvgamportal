## Forecastability across time at the Portal LTREB

This is a study trying to understand how forecastability varies through time at a long-term study of small mammal community dynamics.

## Running existing models

Models are run using the `sliding_window_fits.R` script from the command line:

```sh
Rscript sliding_window_fits.R
```

or from within R:

```r
source("sliding_window_fits.R")
```

## Configuring model runs

Which models, species, test starts, and training window widths are controlled by `config.yml`.

For models the name must match one of the model names in the `fit_model` function in `models.R`.

One or more sets of species to fit is declared in `species_sets`, e.g.:

```yml
  species_sets:
    - name: all

    - name: single
      species: [DM, DO, DS, PB, PP, PF]
      separately: yes
```

Setting `name` to `all` will include all species in the data. Setting `separately` to `yes` (default is `no`) will fit the model separately to each species instead of fitting it jointly.

## Creating and testing new models

### Add a new model

1. Add a new model file to `models/` that describes the model and fits it in a `fit_<model_name>` function. This function should take only three inputs: `data_train`, `data_test`, and `species_list`.
2. Add a new item to the switch statement in the `fit_model` function in `models.R` that calls the new `fit_<model_name>` function.

### Testing the new model

Run the test_model.R script to test the model.

This can be from the command line, e.g.:

```sh
Rscript test_model.R AR_SINGLE_SPECIES --species DM --test-start 240 --train-win-width 60
```

Or from R, e.g.:

```r
source("test_model.R")
fit <- test_model("AR_SINGLE_SPECIES", species = "DM")
```

### Adding the model to full runs

Add the model to `config.yml` so that it is run during the main runs.

## Using the visualization app

There is a shiny app for visualizing forecast scores through time. It can be run from the shell:

```sh
Rscript app.R
```

Or from R:

```r
shiny::runApp()
```
