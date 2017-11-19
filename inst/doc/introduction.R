## ----load-data-----------------------------------------------------------
data(Titanic, package = "datasets")
str(Titanic)

## ----partial-tables------------------------------------------------------
partial_tables <- margin.table(Titanic, c(2,4,1))
partial_tables

## ----marginal-table------------------------------------------------------
marginal_table <- margin.table(Titanic, c(2,4))
marginal_table

## ----prop-table----------------------------------------------------------
prop.table(marginal_table)

## ----marginal-OR---------------------------------------------------------
library(samplesizeCMH)

odds.ratio(marginal_table)

## ----conditional-OR------------------------------------------------------
apply(partial_tables, 3, odds.ratio)

## ----mantelhaen----------------------------------------------------------
library(stats)

mantelhaen.test(partial_tables)

## ----BreslowDay----------------------------------------------------------
library(DescTools)

BreslowDayTest(x = partial_tables, OR = 1)

