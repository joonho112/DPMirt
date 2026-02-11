# ============================================================================
# Global Variable Declarations
# ============================================================================
# Suppress R CMD check NOTEs for NIMBLE DSL variables.
# These variables are used inside nimbleCode() and nimbleFunction() definitions
# and are valid within NIMBLE's domain-specific language (DSL).
# ============================================================================

utils::globalVariables(c(
  # ----- NIMBLE model code variables (used inside nimbleCode blocks) -----
  "M", "N", "delta", "eta",
  "beta.tmp", "gamma.tmp", "lambda", "logLambda.tmp",
  "muTilde", "s2Tilde", "zi",

  # ----- NIMBLE special assignment operators -----
  "logit<-", "log<-",

  # ----- nimbleFunction setup/run variables -----
  "model", "target", "target1", "target2",
  "centering_mean", "mvSaved", "calcNodesNoSelf",
  "adaptive", "adaptiveProcedure", "adaptInterval",
  "timesRan", "timesAccepted", "timesAdapted",
  "gamma1", "optimalAR", "scaleOriginal",
  "nodesToCenter", "numPairs", "samplers", "nodes",
  "copy", "returnType",

  # ----- NIMBLE DSL functions that shadow base R -----
  "dbinom", "rbinom", "rexp",

  # ----- Used in nimbleFunction run code -----
  "x",

  # ----- ggplot2 .data pronoun (used in R/plot_gg.R) -----
  ".data"
))
