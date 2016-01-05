#!env bash

Rscript arff2csv.R -in ../data/ThoraricSurgery.arff -out ../data/ThoraricSurgery.csv

# feature engineering
Rscript feature_engineering.R 

