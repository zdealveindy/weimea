# test_fourth
library (weimea)

# define variables
COM <- vltava$herbs$spe
TRAITS <- vltava$herbs$traits
ENV <- vltava$env$COVERE32
CWM <- cwm (com = COM, traits = TRAITS)

# do some pre-calculation for comparison
COM_p_plant.height <- COM[,!is.na (TRAITS$plant.height)]/rowSums (COM[,!is.na (TRAITS$plant.height)])
COM_p_SLA <- COM[,!is.na (TRAITS$SLA)]/rowSums (COM[,!is.na (TRAITS$SLA)])
COM_p_seed.weight <- COM[,!is.na (TRAITS$seed.weight)]/rowSums (COM[,!is.na (TRAITS$seed.weight)])

calc_cwm_plant.height <- colSums (t (COM_p_plant.height) * (TRAITS$plant.height[!is.na (TRAITS$plant.height)]))
calc_cwm_SLA <- colSums (t (COM_p_SLA) * (TRAITS$SLA[!is.na (TRAITS$SLA)]))
calc_cwm_seed.weight <- colSums (t (COM_p_seed.weight) * (TRAITS$seed.weight[!is.na (TRAITS$seed.weight)]))


test_that ("test_fourth returns correct r_fc",{
  expect_equal (coef (test_fourth (env = ENV, cwm = CWM, perm = 0))$r_fc, c(0.16215121, 0.16271554, -0.02397019))
  expect_equal (coef (test_fourth (env = ENV, cwm = CWM, perm = 0))$r_fc, fourth.corner.ade(sitspe = COM, speatt = TRAITS, env = ENV, fc.test = 6, perm = 0, chessel = FALSE)$fourthcorner)
  expect_equal (coef (test_fourth (env = ENV, cwm = CWM, perm = 0))$r_fc, { # comparison with ade4::fourthcorner
    fc_ade_plant.height <- ade4::fourthcorner (tabR = as.data.frame (ENV), tabL = as.data.frame (COM)[,!is.na (TRAITS$plant.height)], tabQ = as.data.frame (TRAITS$plant.height)[!is.na (TRAITS$plant.height),, drop = F], nrepet = 0)
    fc_ade_SLA <- ade4::fourthcorner (tabR = as.data.frame (ENV), tabL = as.data.frame (COM)[,!is.na (TRAITS$SLA)], tabQ = as.data.frame (TRAITS$SLA)[!is.na (TRAITS$SLA),, drop = F], nrepet = 0)
    fc_ade_seed.weight <- ade4::fourthcorner (tabR = as.data.frame (ENV), tabL = as.data.frame (COM)[,!is.na (TRAITS$seed.weight)], tabQ = as.data.frame (TRAITS$seed.weight)[!is.na (TRAITS$seed.weight),, drop = F], nrepet = 0)
    c(fc_ade_plant.height$tabD$obs, fc_ade_SLA$tabD$obs, fc_ade_seed.weight$tabD$obs)
  })
})

test_that ("cwm is calculated correctly", {
  expect_equal (CWM$plant.height, unname (calc_cwm_plant.height))
  expect_equal (CWM$SLA, unname (calc_cwm_SLA))
  expect_equal (CWM$seed.weight, unname (calc_cwm_seed.weight))
})

