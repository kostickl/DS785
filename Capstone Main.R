############################################################################
#Program
#Purpose: Advanced heart failure cohort
#Programmer: Lee Kostick
#Date: 12.1.2020
############################################################################


impala_Connect()


library(dplyr)
library(tibble)
library(lubridate)
library(ggplot2)
library(stringr)
library(data.table)
library(bit64)
library(bit)
library(fasttime)

#To store tables on the server
source("/user/kostil2/Projects/lee-utils/BuildTmpTable.R")
#Start: functions to pull LVEF values
source("/user/kostil2/Projects/VAD_Devices/lvef_pull.R")
source("/user/kostil2/Projects/VAD_Devices/son_of_clem.R")
source("/user/kostil2/Projects/VAD_Devices/lvef_util.R")
#End: functions to pull LVEF values

#File for lookup codes
source("/user/kostil2/Projects/VAD_Devices/code_lists.R")

#Lookup period for +/- days from NYHA classification
lookup_interval = 45
#Min LVEF value
lvef_threshold = 35
#Days from index date for stoke outcomes
stroke_lookforward = 180

############################################################################
#
#Start: Create look-up tables
#
#NYHA Patient List
#Inotrope Table
#Balloon Pump Table
#HT/WT/BMI Table
#Patient gender/age
#LVEF Function
#Diagnoses
#Temporary LVADS
#Stroke outcomes
#Cardiac flags from updated GLIDE Medtronic table
#
############################################################################


#Get a distinct list of patient ids with a NYHA score of 3 or 4
nyha_base_tbl <-
  tbl(im.con, 'cvg_device_ehr_2019.tbl_cde_hf_nyha') %>%
  filter(nyha_score >= 3) %>%
  mutate(nyha_date = to_date(nyha_date)) %>%
  group_by(patient_id) %>%
  summarise(min_nyha_dt = min(nyha_date))

nyha_pat_list = nyha_base_tbl %>% select(patient_id)

#Pull patients with at least one inotrope administered
inotrope_table <-
  tbl(im.con, 'cvg_device_ehr_2019.ehr_rx_adm') %>%  
  filter(generic_desc %in% inotrope_names & !is.na(admin_date_time)) %>%
  mutate(admin_date = to_date(admin_date_time)) %>%
  select(patient_id, admin_date, generic_desc) %>%
  inner_join(nyha_pat_list, by="patient_id")

#Pull patients with at least one balloon pump proc code
balloon_table <-
  tbl(im.con, 'cvg_device_ehr_2019.ehr_proc') %>%  
  filter(proc_code %in% iabp_codes) %>%
  select(patient_id, proc_date) %>%
  inner_join(nyha_pat_list, by="patient_id")


#Grab ht, wt, and bmi values
ht_wt_table_pre <-
  tbl(im.con, 'cvg_device_ehr_2019.ehr_obs') %>%  
  filter(obs_type %in% c('HT','WT', 'BMI')) %>%
  select(patient_id, obs_type, obs_unit, obs_result, obs_date) %>%
  inner_join(nyha_pat_list, by='patient_id')

#Reguire heights to be within 5% of the patient's median height, wieght within 40%
#This eliminates data entry issues
ht_wt_table =
  ht_wt_table_pre %>%
  group_by(patient_id, obs_type) %>%
  summarise(obs_pat_median = appx_median(obs_result, na.rm=TRUE)) %>%
  inner_join(ht_wt_table_pre, by=c('patient_id', 'obs_type')) %>%
  mutate(obs_pat_median = as.numeric(obs_pat_median)) %>%
  filter((obs_result == 'HT' & between(as.numeric(obs_result),.95*obs_pat_median, 1.05*obs_pat_median)) |
           (between(as.numeric(obs_result),.6*obs_pat_median, 1.4*obs_pat_median))) %>%
  BuildTmpTable(im.con, 'ht_wt_table', 'kostickl')

#Pull patient gender and birth yr
pat_table <-
  tbl(im.con, 'cvg_device_ehr_2019.ehr_patient') %>%  
  select(patient_id, birth_yr, gender) %>%
  mutate(birth_yr = as.numeric(substr(birth_yr,1,4))) %>%
  inner_join(nyha_pat_list, by='patient_id') %>%
  BuildTmpTable(im.con, 'pat_table', 'kostickl')


#To pass to LVEF function
pat_list_collect = nyha_base_tbl %>% collect()

#Internal function to pull LVEF values
lvef_table = lvef_faust(cohort = pat_list_collect)

lvef_table_filt <- 
  lvef_table %>% 
  filter(raw_val < lvef_threshold)

#Pull diagnoses 
dx_table <-
  tbl(im.con, 'cvg_device_ehr_2019.ehr_diag') %>%
  filter(diagnosis_cd_type %in% c('ICD10','ICD9')) %>%
  filter(diagnosis_cd %in% hf_codes |
           diagnosis_cd %in% cm_codes |
           diagnosis_cd %in% chd_codes |
           diagnosis_cd %in% af_codes |
           diagnosis_cd %in% stroke_codes |
           substr(diagnosis_cd,1,3) %in% ami_codes) %>%
  inner_join(nyha_pat_list, by='patient_id') %>%
  select(patient_id, diag_date_time, diagnosis_cd, diag_date_time) %>%
  mutate(diag_dt = to_date(diag_date_time),
         dx_cat = case_when(diagnosis_cd %in% hf_codes |
                              diagnosis_cd %in% cm_codes ~ 'HF',
                            diagnosis_cd %in% chd_codes ~ 'CHD',
                            diagnosis_cd %in% af_codes ~ 'AF',
                            diagnosis_cd %in% stroke_codes ~ 'Stroke',
                            substr(diagnosis_cd,1,3) %in% ami_codes ~ 'AMI')) %>%
  select(-diag_date_time, -diagnosis_cd) %>%
  BuildTmpTable(im.con, 'dx_table', 'kostickl')

#Temporary LVADs
temp_lvads = 
  tbl(im.con, 'cvg_device_ehr_2019.ehr_proc') %>%  
  filter(proc_code %in% temp_lvad_codes) %>%
  select(patient_id, proc_date) %>%
  inner_join(nyha_pat_list, by="patient_id") %>%
  group_by(patient_id) %>%
  summarise(first_lvad = min(proc_date))

#Stroke outcomes
outcomes <- 
  tbl(im.con, 'cvg_device_ehr_2019.ehr_diag') %>%  
  filter(diagnosis_cd %in% ischemic_stroke_codes |
         diagnosis_cd %in% hemor_stroke_codes) %>%
  inner_join(nyha_base_tbl, by='patient_id') %>%
  select(diag_date, patient_id) %>%
  collect()

#Covariates from updated GLIDE table; internal to Medtronic
ugli_tab <- 
  tbl(im.con, 'cvg_device_ehr_2019.tbl_hf_ugli') %>%  
  select(patient_id,
         race, 
         region, 
         death, 
         atrial_fibrillation, 
         angina,
         av_block_highdeg, 
         cardmyo_hypertroph, 
         cardmyo_nonischemic, 
         hypertension, 
         sick_sinus_syndrome,
         vt_any, 
         syncope,
         dev_tavr,
         ablation_av_node, 
         ablation_af_pvi, 
         ablation_vt,
         chronic_kidney_disease,
         dialysis, 
         peripheral_vasc, 
         diabetes, 
         apnea,
         copd,
         depression,
         dementia,
         rx_betablocker_any, 
         rx_ace_arb,
         rx_mra,
         rx_diur_loop, 
         rx_diur_thiazide, 
         rx_entresto,
         rx_ap_aspirin,
         rx_ap_antiplatelets,
         rx_ac_noac,
         rx_ac_warfarin,
         rx_anti_arrhythmic_class_i, 
         rx_anti_arrhythmic_amiodarone, 
         rx_sglt_2,
         rx_ca_channel_blocker,
         sds_cardmyo_dilated, 
         sds_collapse,
         sds_dyssynchrony,
         sds_ischemic,
         sds_left_bundle,
         sds_vt_other,
         vt_fibrillation,
         vt_sustained,
         coronary_pci,
         aicd_ind,
         crtd_ind_block_hf) %>% 
  inner_join(nyha_base_tbl, by='patient_id') %>%
  collect()


ugli_tab %>% view()
############################################################################
#
#End: Create look-up tables
#
############################################################################


#############################################
#Build table combining
# 1. inotropes
# 2. balloon pumps
# 3. ht, wt, BMI
# 4. Previous temporary lvad
# 5. Diagnosis flags
# 6. LVEF values
#############################################

#Table with all NYHA 3 or 4 patients
#Include first inotrope use that is within lookup_inerval
combined1 <-
  nyha_base_tbl %>%
  #Add birth_yr and gender
  left_join(pat_table, by="patient_id") %>%
  mutate(age = year(min_nyha_dt) - birth_yr) %>%
  #Add inotropes
  left_join(inotrope_table, by="patient_id") %>%
  mutate(days_diff = abs(datediff(admin_date, min_nyha_dt))) %>%
  arrange(patient_id, days_diff) %>%
  group_by(patient_id, age) %>%
  mutate(first_ino = min(if_else(days_diff <= lookup_interval, admin_date, '9999-12-31')),
         ino_cnt = n()) %>%
  filter(row_number()==1) %>%
  ungroup() 

#Insert dates for balloon pump
combined2 <-
  combined1 %>%
  #Add balloon pump
  left_join(balloon_table, by="patient_id") %>%
  mutate(days_diff = abs(datediff(proc_date, min_nyha_dt))) %>%
  arrange(patient_id, days_diff) %>%
  group_by(patient_id) %>%
  mutate(first_balloon = min(if_else(days_diff <= lookup_interval, proc_date, '9999-12-31'))) %>%
  filter(row_number()==1) %>%
  ungroup()

#Add ht, wt, and bmi
combined3 <-
  combined2 %>%
  left_join(ht_wt_table, by = "patient_id") %>%
  mutate(days_diff = abs(datediff(obs_date, min_nyha_dt))) %>%
  arrange(patient_id, obs_type, days_diff) %>%
  group_by(patient_id, obs_type) %>%
  filter(row_number()==1) %>%
  ungroup() %>%
  mutate(ht = if_else(obs_type == "HT", obs_result, NULL),
         ht_dt = if_else(obs_type == "HT", to_date(obs_date), NULL),
         wt = if_else(obs_type == "WT", obs_result, NULL),
         wt_dt = if_else(obs_type == "WT", to_date(obs_date), NULL),
         bmi = if_else(obs_type == "BMI", obs_result, NULL),
         bmi_dt = if_else(obs_type == "BMI", to_date(obs_date), NULL)) %>%
  group_by(patient_id, min_nyha_dt, first_ino, ino_cnt, first_balloon, birth_yr, gender, age) %>%
  summarise(ht = max(ht), ht_dt = max(ht_dt),
            wt = max(wt), wt_dt = max(wt_dt),
            bmi = max(bmi), bmi_dt = max(bmi_dt)) 

#Temporary lvads
combined4 <-
  combined3 %>%
  left_join(temp_lvads, by="patient_id") %>%
  mutate(lvad_is_before = if_else(first_lvad < min_nyha_dt & !is.na(first_lvad), 1,0))

#Add dx flags
combined5 <-
  combined4 %>%
  left_join(dx_table, by='patient_id') %>%
  mutate(days_diff = abs(datediff(diag_dt, min_nyha_dt))) %>%
  arrange(patient_id, dx_cat, days_diff) %>%
  group_by(patient_id, dx_cat) %>%
  filter(row_number()==1) %>%
  ungroup() %>%
  mutate(hf_dt = if_else(dx_cat == "HF", diag_dt, NULL),
         chd_dt = if_else(dx_cat == "CHD", diag_dt, NULL),
         af_dt = if_else(dx_cat == "AF", diag_dt, NULL),
         stroke_dt = if_else(dx_cat == "Stroke", diag_dt, NULL),
         ami_dt = if_else(dx_cat == "AMI", diag_dt, NULL)) %>%
  group_by(patient_id, min_nyha_dt, first_ino, ino_cnt, first_balloon,
           ht, wt, bmi, ht_dt, wt_dt, bmi_dt, birth_yr, gender, age, lvad_is_before) %>%
  summarise(hf_dt = max(hf_dt), chd_dt = max(chd_dt), af_dt = max(af_dt),
            stroke_dt = max(stroke_dt), ami_dt = max(ami_dt)) %>%
  collect()

#Add LVEF values
combined6 <-
  combined5 %>%
  left_join(lvef_table_filt, by="patient_id") %>%
  mutate(days_diff = abs(difftime(EVENT_DATE, min_nyha_dt, units = "days"))) %>%
  arrange(patient_id, days_diff) %>%
  group_by(patient_id) %>%
  filter(row_number()==1) %>%
  mutate(first_lvef = case_when(days_diff <= lookup_interval ~ EVENT_DATE, 
                                is.na(days_diff) ~ as.Date('9999-12-30'),
                                TRUE ~ as.Date('9999-12-31'))) %>%
  
  ungroup() %>%
  select(-EVENT_DATE, -raw_val, -mod, -days_diff, -EVENT_NAME)

#BSA calculation
combined6$bsa = sqrt(((as.numeric(combined6$ht)*as.numeric(combined6$wt))/3600))

#Apply inclusion/exclusion criteria
base_inc = 
  combined6 %>%
  filter(min_nyha_dt >= '2015-01-01') %>%
  mutate(include = if_else(
         age >= 18 & !is.na(age) &
         bsa > 1.2 & !is.na(bsa) &
         lvad_is_before == 0 &
         ((!is.na(first_ino) & year(first_ino) != 9999) | (!is.na(first_balloon) & year(first_balloon) != 9999)) &
         year(first_lvef) != 9999 &
         bmi < 40 & !is.na(bmi) &
         is.na(chd_dt) &
         (is.na(ami_dt) | !between(difftime(ami_dt, min_nyha_dt, units = 'days'), -30, 0)),
         1,0))

#Add index date, first occurrence of lvef, inotrope, iabl, or NYHA
base_ind =
  base_inc %>% 
  filter(include == 1) %>%
  group_by(patient_id) %>%
  mutate(index_dt = min(min_nyha_dt,
                        if_else(!is.na(first_ino), first_ino, '9999-12-31'),
                        if_else(!is.na(first_balloon), as.character(first_balloon), '9999-12-31'),
                        if_else(!is.na(first_lvef), as.character(first_lvef), '9999-12-31')))

base_ind %>% 
  filter(!is.na(first_balloon) & is.na(first_ino)) %>%
  group_by() %>%
  summarise(n_distinct(patient_id))

#Add indicator for stroke within 
#4243 patients
base_str =
  base_ind %>%
  left_join(outcomes, by='patient_id') %>%
  group_by(patient_id) %>%
  mutate(had_stroke = if_else(between(difftime(diag_date, index_dt, units = 'days'), 0, stroke_lookforward) & !is.na(diag_date),1,0),
         old_stroke = if_else(between(difftime(diag_date, index_dt, units = 'days'), -3000, -90) & !is.na(diag_date),1,0)) %>%
  group_by(patient_id,
           ino_cnt,
           bmi,
           gender,
           age,
           lvad_is_before,
           stroke_dt,
           bsa,
           index_dt) %>%
  summarise(had_stroke = max(had_stroke), old_stroke = max(old_stroke)) %>%
  ungroup()

#Inner join to UGLI table for list of covariates
#3926 patients
base_str_cov =
  base_str %>%
  inner_join(ugli_tab, by='patient_id') %>% 
  mutate(atrial_fibrillation_flag =  if_else(atrial_fibrillation < index_dt & (!is.na(atrial_fibrillation) &  !is.null(atrial_fibrillation)), 1,0),
                                             angina_flag =  if_else(angina < index_dt & (!is.na(angina) &  !is.null(angina)), 1,0),
                                             av_block_highdeg_flag =  if_else(av_block_highdeg < index_dt & (!is.na(av_block_highdeg) &  !is.null(av_block_highdeg)), 1,0),
                                             cardmyo_hypertroph_flag =  if_else(cardmyo_hypertroph < index_dt & (!is.na(cardmyo_hypertroph) &  !is.null(cardmyo_hypertroph)), 1,0),
                                             cardmyo_nonischemic_flag =  if_else(cardmyo_nonischemic < index_dt & (!is.na(cardmyo_nonischemic) &  !is.null(cardmyo_nonischemic)), 1,0),
                                             hypertension_flag =  if_else(hypertension < index_dt & (!is.na(hypertension) &  !is.null(hypertension)), 1,0),
                                             sick_sinus_syndrome_flag =  if_else(sick_sinus_syndrome < index_dt & (!is.na(sick_sinus_syndrome) &  !is.null(sick_sinus_syndrome)), 1,0),
                                             vt_any_flag =  if_else(vt_any < index_dt & (!is.na(vt_any) &  !is.null(vt_any)), 1,0),
                                             syncope_flag =  if_else(syncope < index_dt & (!is.na(syncope) &  !is.null(syncope)), 1,0),
                                             dev_tavr_flag =  if_else(dev_tavr < index_dt & (!is.na(dev_tavr) &  !is.null(dev_tavr)), 1,0),
                                             ablation_av_node_flag =  if_else(ablation_av_node < index_dt & (!is.na(ablation_av_node) &  !is.null(ablation_av_node)), 1,0),
                                             ablation_af_pvi_flag =  if_else(ablation_af_pvi < index_dt & (!is.na(ablation_af_pvi) &  !is.null(ablation_af_pvi)), 1,0),
                                             ablation_vt_flag =  if_else(ablation_vt < index_dt & (!is.na(ablation_vt) &  !is.null(ablation_vt)), 1,0),
                                             chronic_kidney_disease_flag =  if_else(chronic_kidney_disease < index_dt & (!is.na(chronic_kidney_disease) &  !is.null(chronic_kidney_disease)), 1,0),
                                             dialysis_flag =  if_else(dialysis < index_dt & (!is.na(dialysis) &  !is.null(dialysis)), 1,0),
                                             peripheral_vasc_flag =  if_else(peripheral_vasc < index_dt & (!is.na(peripheral_vasc) &  !is.null(peripheral_vasc)), 1,0),
                                             diabetes_flag =  if_else(diabetes < index_dt & (!is.na(diabetes) &  !is.null(diabetes)), 1,0),
                                             apnea_flag =  if_else(apnea < index_dt & (!is.na(apnea) &  !is.null(apnea)), 1,0),
                                             copd_flag =  if_else(copd < index_dt & (!is.na(copd) &  !is.null(copd)), 1,0),
                                             depression_flag =  if_else(depression < index_dt & (!is.na(depression) &  !is.null(depression)), 1,0),
                                             dementia_flag =  if_else(dementia < index_dt & (!is.na(dementia) &  !is.null(dementia)), 1,0),
                                             rx_betablocker_any_flag =  if_else(rx_betablocker_any < index_dt & (!is.na(rx_betablocker_any) &  !is.null(rx_betablocker_any)), 1,0),
                                             rx_ace_arb_flag =  if_else(rx_ace_arb < index_dt & (!is.na(rx_ace_arb) &  !is.null(rx_ace_arb)), 1,0),
                                             rx_mra_flag =  if_else(rx_mra < index_dt & (!is.na(rx_mra) &  !is.null(rx_mra)), 1,0),
                                             rx_diur_loop_flag =  if_else(rx_diur_loop < index_dt & (!is.na(rx_diur_loop) &  !is.null(rx_diur_loop)), 1,0),
                                             rx_diur_thiazide_flag =  if_else(rx_diur_thiazide < index_dt & (!is.na(rx_diur_thiazide) &  !is.null(rx_diur_thiazide)), 1,0),
                                             rx_entresto_flag =  if_else(rx_entresto < index_dt & (!is.na(rx_entresto) &  !is.null(rx_entresto)), 1,0),
                                             rx_ap_aspirin_flag =  if_else(rx_ap_aspirin < index_dt & (!is.na(rx_ap_aspirin) &  !is.null(rx_ap_aspirin)), 1,0),
                                             rx_ap_antiplatelets_flag =  if_else(rx_ap_antiplatelets < index_dt & (!is.na(rx_ap_antiplatelets) &  !is.null(rx_ap_antiplatelets)), 1,0),
                                             rx_ac_noac_flag =  if_else(rx_ac_noac < index_dt & (!is.na(rx_ac_noac) &  !is.null(rx_ac_noac)), 1,0),
                                             rx_ac_warfarin_flag =  if_else(rx_ac_warfarin < index_dt & (!is.na(rx_ac_warfarin) &  !is.null(rx_ac_warfarin)), 1,0),
                                             rx_anti_arrhythmic_class_i_flag =  if_else(rx_anti_arrhythmic_class_i < index_dt & (!is.na(rx_anti_arrhythmic_class_i) &  !is.null(rx_anti_arrhythmic_class_i)), 1,0),
                                             rx_anti_arrhythmic_amiodarone_flag =  if_else(rx_anti_arrhythmic_amiodarone < index_dt & (!is.na(rx_anti_arrhythmic_amiodarone) &  !is.null(rx_anti_arrhythmic_amiodarone)), 1,0),
                                             rx_sglt_2_flag =  if_else(rx_sglt_2 < index_dt & (!is.na(rx_sglt_2) &  !is.null(rx_sglt_2)), 1,0),
                                             rx_ca_channel_blocker_flag =  if_else(rx_ca_channel_blocker < index_dt & (!is.na(rx_ca_channel_blocker) &  !is.null(rx_ca_channel_blocker)), 1,0),
                                             sds_cardmyo_dilated_flag =  if_else(sds_cardmyo_dilated < index_dt & (!is.na(sds_cardmyo_dilated) &  !is.null(sds_cardmyo_dilated)), 1,0),
                                             sds_collapse_flag =  if_else(sds_collapse < index_dt & (!is.na(sds_collapse) &  !is.null(sds_collapse)), 1,0),
                                             sds_dyssynchrony_flag =  if_else(sds_dyssynchrony < index_dt & (!is.na(sds_dyssynchrony) &  !is.null(sds_dyssynchrony)), 1,0),
                                             sds_ischemic_flag =  if_else(sds_ischemic < index_dt & (!is.na(sds_ischemic) &  !is.null(sds_ischemic)), 1,0),
                                             sds_left_bundle_flag =  if_else(sds_left_bundle < index_dt & (!is.na(sds_left_bundle) &  !is.null(sds_left_bundle)), 1,0),
                                             sds_vt_other_flag =  if_else(sds_vt_other < index_dt & (!is.na(sds_vt_other) &  !is.null(sds_vt_other)), 1,0),
                                             vt_fibrillation_flag =  if_else(vt_fibrillation < index_dt & (!is.na(vt_fibrillation) &  !is.null(vt_fibrillation)), 1,0),
                                             vt_sustained_flag =  if_else(vt_sustained < index_dt & (!is.na(vt_sustained) &  !is.null(vt_sustained)), 1,0),
                                             coronary_pci_flag =  if_else(coronary_pci < index_dt & (!is.na(coronary_pci) &  !is.null(coronary_pci)), 1,0),
                                             aicd_ind_flag =  if_else(aicd_ind < index_dt & (!is.na(aicd_ind) &  !is.null(aicd_ind)), 1,0),
                                             crtd_ind_block_hf_flag =  if_else(crtd_ind_block_hf < index_dt & (!is.na(crtd_ind_block_hf) &  !is.null(crtd_ind_block_hf)), 1,0)) %>%
  #Add flag for death outcome
  mutate(had_death = if_else(between(difftime(death, index_dt, units = 'days'), 0, 180) & !is.na(death),1,0)) %>%
  group_by(patient_id) %>%
  mutate(had_strk_or_dth = max(had_stroke, had_death),
         cardmyo_flag = max(cardmyo_hypertroph_flag, cardmyo_nonischemic_flag, rx_anti_arrhythmic_class_i_flag, rx_anti_arrhythmic_amiodarone_flag),
         ablation_flag = max(ablation_av_node_flag, ablation_af_pvi_flag, ablation_vt_flag),
         diurr_flag = max(rx_diur_loop_flag, rx_diur_thiazide_flag),
         dep_dem_flag = max(depression_flag, dementia_flag),
         rx_ap_antiplat_flag = max(rx_ap_aspirin_flag, rx_ap_antiplatelets_flag),
         rx_anticoag_flag = max(rx_ac_noac_flag, rx_ac_warfarin_flag, aicd_ind_flag),
         af_flag = max(atrial_fibrillation_flag, sds_left_bundle_flag),
         ssskidney_flag = max(chronic_kidney_disease_flag, sick_sinus_syndrome_flag),
         diab_cond_flag = max(diabetes_flag, peripheral_vasc_flag), 
         beta_flag = max(rx_ace_arb_flag, rx_betablocker_any_flag),
         vt_flag = max(vt_fibrillation_flag, vt_sustained_flag, vt_any_flag, sds_vt_other_flag)) %>%
  dplyr::select(-atrial_fibrillation,
                -atrial_fibrillation_flag,
                 -angina,
                 -av_block_highdeg,
                 -cardmyo_hypertroph,
                 -cardmyo_nonischemic,
                 -cardmyo_hypertroph_flag, 
                 -cardmyo_nonischemic_flag,
                 -hypertension,
                 -sick_sinus_syndrome,
                 -sick_sinus_syndrome_flag,
                 -vt_any,
                 -vt_any_flag,
                 -syncope,
                 -dev_tavr,
                 -ablation_av_node,
                 -ablation_af_pvi,
                 -ablation_vt,
                 -ablation_av_node_flag,
                 -ablation_af_pvi_flag,
                 -ablation_vt_flag,
                 -chronic_kidney_disease,
                 -chronic_kidney_disease_flag,
                 -dialysis,
                 -peripheral_vasc,
                 -peripheral_vasc_flag,
                 -diabetes,
                 -diabetes_flag,
                 -apnea,
                 -copd,
                 -depression,
                 -dementia,
                 -depression_flag,
                 -dementia_flag,
                 -rx_betablocker_any,
                 -rx_betablocker_any_flag,
                 -rx_ace_arb,
                 -rx_ace_arb_flag,
                 -rx_mra,
                 -rx_diur_loop,
                 -rx_diur_thiazide,
                 -rx_diur_loop_flag,
                 -rx_diur_thiazide_flag,
                 -rx_entresto,
                 -rx_ap_aspirin,
                 -rx_ap_antiplatelets,
                 -rx_ac_noac,
                 -rx_ac_warfarin,
                 -rx_ap_aspirin_flag,
                 -rx_ap_antiplatelets_flag,
                 -rx_ac_noac_flag,
                 -rx_ac_warfarin_flag,
                 -rx_anti_arrhythmic_class_i,
                 -rx_anti_arrhythmic_amiodarone,
                 -rx_anti_arrhythmic_class_i_flag,
                 -rx_anti_arrhythmic_amiodarone_flag,
                 -rx_sglt_2,
                 -rx_ca_channel_blocker,
                 -sds_cardmyo_dilated,
                 -sds_collapse,
                 -sds_dyssynchrony,
                 -sds_ischemic,
                 -sds_left_bundle,
                 -sds_left_bundle_flag,
                 -sds_vt_other,
                 -sds_vt_other_flag,
                 -vt_fibrillation,
                 -vt_sustained,
                 -vt_fibrillation_flag,
                 -vt_sustained_flag,
                 -coronary_pci,
                 -aicd_ind,
                 -aicd_ind_flag,
                 -crtd_ind_block_hf,
                 -lvad_is_before)


###############################################################
#
#Check for balanced covariates between cohorts (stoke vs non stroke)
#
###############################################################

data.frame(colnames(base_str_cov))

tot_cnt = nrow(base_str_cov)
tot_cnt_str = nrow(base_str_cov[which(base_str_cov$had_strk_or_dth==1),])
tot_cnt_nostr = nrow(base_str_cov[which(base_str_cov$had_strk_or_dth==0),])
tot_cnt_str
cols = c(4,seq(15,46))

for(i in cols){
  
  print(colnames(base_str_cov[,i]))
  
  print(table(base_str_cov[,i])/tot_cnt)
  
  print("Had stroke or death")
  print(table(base_str_cov[which(base_str_cov$had_strk_or_dth==1),i])/tot_cnt_str)
  
  print("No stroke or death")
  print(table(base_str_cov[which(base_str_cov$had_strk_or_dth==0),i])/tot_cnt_nostr)
}

###############################################################
#
#End - Check for balanced covariates between cohorts (stoke vs non stroke)
#
###############################################################

final_base_table =
  base_str_cov %>%
  dplyr::select(-av_block_highdeg_flag, -dev_tavr_flag, -dialysis_flag, -copd_flag, -rx_entresto_flag, -rx_sglt_2_flag, -sds_collapse_flag,
         -sds_dyssynchrony_flag, -crtd_ind_block_hf_flag,
         -ablation_flag)
