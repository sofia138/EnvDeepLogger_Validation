# =======================================================================
# Analysing performance of temperature low-cost loggers: 2 types and 3 units each

#EnvDeep01 (*6515 0F, tempMad01, tempMad02, tempMad04)
#EnvDeep02 (*F202 08, tempMad03)
#EnvDeep03 (*E31E 0F, tempMad05)
#EnvLogger01 (7A1D 08, tempMad06)
#EnvLogger02 (*F622 0F, tempMad07)
#EnvLogger03 (*C519 0E, tempMad08)
# =======================================================================

library(dplyr) #as.data.frame(), group_by(), mutate(), filter(), arrange(), summarise()
library(tidyr) #crossing
library(ggplot2)
library(lubridate) #time parsing (ymd_hms)
library(data.table) #roll = "nearest"
library(zoo) #rollapply, rollmean
library(dunn.test) #Kruskal-Wallis tests
library(stringr) #str_detect
library(patchwork) #plot_layout

setwd("C:/Users/asnog/Documents/PhD-MareMadeira/TemperatureLogger/")

#. =========================================================================
# 1. Loading loggers' data -----------------------------------------------------
#. =========================================================================

read_logger <- function(file, file_id, logger) {
  read.csv(file, sep = ",", header = TRUE) %>%
    filter(!is.na(temp)) %>%
    mutate(
      time = lubridate::ymd_hms(time, tz = "UTC"),
      file_id = file_id,
      logger = logger
    )
}

logger_files <- tibble::tibble(
  file = c(
    "TemperaturesExpeditionDeepSea2024/tempMad01.csv",
    "TemperaturesExpeditionDeepSea2024/tempMad02.csv",
    "TemperaturesExpeditionDeepSea2024/tempMad03.csv",
    "Scripts&Outputs/EnvDeep01_03_T7.3-01-03_vs_CTD_oom/logger_files/tempMad04.csv",
    "Scripts&Outputs/EnvDeep01_03_T7.3-01-03_vs_CTD_oom/logger_files/tempMad05.csv",
    "Scripts&Outputs/EnvDeep01_03_T7.3-01-03_vs_CTD_oom/logger_files/tempMad06.csv",
    "Scripts&Outputs/EnvDeep01_03_T7.3-01-03_vs_CTD_oom/logger_files/tempMad07.csv",
    "Scripts&Outputs/EnvDeep01_03_T7.3-01-03_vs_CTD_oom/logger_files/tempMad08.csv"
  ),
  file_id = c( # keeps unique ID of logger files aka dives
    "tempMad01", "tempMad02", "tempMad03",
    "tempMad04", "tempMad05", 
    "tempMad06", "tempMad07", "tempMad08"
  ),
  logger = c( # logger ids
    "EnvDeep01", "EnvDeep01", "EnvDeep02",
    "EnvDeep01", "EnvDeep03",
    "EnvLogger01", "EnvLogger02", "EnvLogger03"
  )
)

loggers_all <- purrr::pmap_dfr(
  logger_files,
  read_logger
)

loggers_all <- loggers_all %>%
  select(dive = file_id, time, log_temp=temp, logger)

#.------------------------------------------------------------------------------
# Exploratory - Visualization of the temperature per time ----------------------
#.------------------------------------------------------------------------------
# Plot function
plot_logger_ts <- function(df, dive_name) {
  df %>% 
    filter(dive == dive_name) %>%
    ggplot(aes(x = time, y = log_temp)) +
    geom_line(linewidth = 0.6, color = "dodgerblue3") +
    theme_minimal(base_size = 14) +
    scale_x_datetime(
      date_breaks = "6 hours",
      date_minor_breaks = "1 hour",
      date_labels = "%b %d\n%H:%M"
    ) +
    labs(
      title = paste("Temperature Time Series –", dive_name),
      x = "Time (UTC)",
      y = "Temperature (°C)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.6),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.4),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Apply to all
unique(loggers_all$dive) %>%
  lapply(function(x) plot_logger_ts(loggers_all, x))

#. ==============================================================================
# 2. Loading and Pre-processing data from reference instruments -----------------
#. ==============================================================================

#.-------------------------------------------------------------------------------------------
# 2.1. Preparing functions for cleaning near-surface instability for D01-D03 ----------------
#.-------------------------------------------------------------------------------------------

## depth_threshold - min depth to count for the dive
## window_n - how many future points must stay below threshold

#clean initial instability
clean_start_depth <- function(df, depth_col = "p", time_col = "time",
                              depth_threshold = 10, window_n = 100) {
  df <- df %>%
    arrange(.data[[time_col]]) %>%
    mutate(
      dp = c(diff(.data[[depth_col]]), NA) # depth change
    )
  df <- df %>%
    mutate(
      # condition 1: future always below threshold
      future_min = zoo::rollapply(
        .data[[depth_col]],
        width = window_n,
        FUN = min,
        fill = NA,
        align = "left"
      ),
      # condition 2: mostly increasing (allow small noise)
      future_increasing = zoo::rollapply(
        dp > -0.5,   # allow small negative wiggles
        width = window_n,
        FUN = function(x) mean(x, na.rm = TRUE) > 0.9,
        fill = NA,
        align = "left"
      )
    )
  # find first valid time
  t_start_clean <- df %>%
    filter(
      .data[[depth_col]] > depth_threshold,
      future_min > depth_threshold,
      future_increasing == TRUE
    ) %>%
    slice(1) %>%
    pull(.data[[time_col]])
  # subset
  df %>%
    filter(.data[[time_col]] >= t_start_clean) %>%
    select(-dp, -future_min, -future_increasing)
}

#Clean ascent - remove anything shallower than 10 m
clean_end_ascent <- function(df,
                             depth_col = "p",
                             time_col = "time",
                             depth_threshold = 10) { df <- df %>%
                               arrange(.data[[time_col]])
                             t_cut <- df %>% # find the earliest time when its deeper than 10 m
                               filter(.data[[depth_col]] >= depth_threshold) %>%
                               summarise(last_time = max(.data[[time_col]], na.rm = TRUE)) %>%
                               pull(last_time)
                             if (is.infinite(t_cut) || is.na(t_cut)) { # safety check but should not happen
                               warning("No depth >= threshold found — returning full dataset")
                               return(df)
                             } # trim everything after that point
                             df %>%
                               filter(.data[[time_col]] <= t_cut)
}

#.------------------------------------------------------------------------------------
# 2.2. Read and Process ROV file for D01 ---------------------------------------------------------
#.------------------------------------------------------------------------------------

# Read file
read_sbe49 <- function(filepath) {
  raw <- readLines(filepath) #Read entire file
  data_start <- which(grepl("^\\s*\\d+\\s+\\d", raw))[1] # Find first line that looks like numeric data
  header <- raw[1:data_start] #Read the header for names
  names_line <- header[grepl("^# name", header)]
  colnames <- sub("^# name \\d+ =\\s*", "", names_line) # Extract column names from "# name X = variable"
  idx_time  <- which(colnames == "timeJ: Julian Days")
  idx_depth <- which(colnames == "prdM: Pressure, Strain Gauge [db]")
  idx_temp  <- which(colnames == "t4990C: Temperature [ITS-90, deg C]")
  needed_idx <- c(idx_time, idx_depth, idx_temp)
  df <- read.table(text = raw[data_start:length(raw)], header = FALSE)[, needed_idx]
  names(df) <- c("time", "p", "t")
  df
}

# Convert date
convert_julianday <- function(jd) {
  origin <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC") #assumes this date as the origin
  origin + (jd - 1) * 86400
}

# Process ctd data and return a data frame
process_rov <- function(filepath) {
  df <- read_sbe49(filepath)
  df |>
    mutate(
      time = convert_julianday(time), #converts Julian days (numeric) to POSIXct datetime
      p = p, #assumes pressure (db) ~ m
      temp  = t #temperature (°C)
    ) |>
    select(time, p, temp)|>
    # Cleaning Steps
    clean_start_depth(
      depth_col = "p", 
      time_col = "time", 
    ) |>
    clean_end_ascent(
      depth_col = "p", 
      time_col = "time", 
    )
}

df_tmp <- read_sbe49("CTDdataMSM126andM209/msm_126_1/rov/msm_126_rov_012_1sec.cnv")
names(df_tmp) #confirms names of the header

rov01 <- process_rov("CTDdataMSM126andM209/msm_126_1/rov/msm_126_rov_012_1sec.cnv")


# Define ROV phases of the dive ----------------------------------------

#Find end of Descent/Initial bottom time
rov_work <- rov01 %>% arrange(time) #sanity check -  ensures data is in chronological order

t_start <- min(rov_work$time)
t_end   <- max(rov_work$time)

dive_duration <- as.numeric(difftime(t_end, t_start, units = "secs"))

#End of descent
idx_descent_end <- which.max(rov_work$p)
t_descent_end <- rov_work$time[idx_descent_end]

#FIND END OF BOTTOM TIME/INITIAL ASCENT
rov_work <- rov_work %>%
  mutate(
    depth_smooth = rollmean(p, k = 61, fill = NA, align = "center")
  )
rov_work <- rov_work %>%
  mutate(
    w = c(NA, diff(depth_smooth)) / c(NA, diff(as.numeric(time)))
  )
rov_rev <- rov_work %>%
  arrange(desc(time))
w_up_thresh <- quantile(rov_rev$w, 0.25, na.rm = TRUE)
rov_rev <- rov_rev %>%
  mutate(
    is_ascending = w < w_up_thresh
  )
dt <- median(diff(as.numeric(rov_rev$time)), na.rm = TRUE)
min_duration <- 300  # seconds
min_n <- ceiling(min_duration / dt)
r <- rle(rov_rev$is_ascending)
ends <- cumsum(r$lengths)
starts <- ends - r$lengths + 1

runs <- data.frame(
  is_ascending = r$values,
  start = starts,
  end   = ends,
  n     = r$lengths
)

asc_run <- runs %>%
  filter(is_ascending, n >= min_n) %>%
  slice(1)
t_ascent_start <- rov_rev$time[asc_run$end]


#Attribute phases to each reading
rov_work <- rov_work %>%
  mutate(
    phase = case_when(
      time >= t_start & time < t_descent_end ~ "Descent",
      time >= t_descent_end & time < t_ascent_start ~ "Bottom",
      time >= t_ascent_start & time <= t_end ~ "Ascent",
      TRUE ~ NA_character_
    )
  )

rov01 <- rov_work %>% mutate(
  dive = "ROV Dive") %>%
  select(dive, time, p, ctd_temp = temp, phase)

#.---------------------------------------------------------------------------------------
# 2.3. Read CTD data for D02 and D03 ----------------------------------------------------
#.---------------------------------------------------------------------------------------

# Matlab datenum → POSIXct
matlab2posix <- function(dn) {
  as.POSIXct((dn - 719529) * 86400, origin = "1970-01-01", tz = "UTC")
}

#Create function to extract tim (converted to time) and pressure from the ctd files
read_ctd_tim_p <- function(filepath) {
  lines <- readLines(filepath)
  first_data_line <- grep("^\\s*\\d", lines)[1] # Find the first line that starts with a number
  df <- read.table(text = lines[first_data_line:length(lines)],
                   header = FALSE,
                   fill = TRUE)
  df <- df[, 1:3]
  names(df) <- c("tim", "p", "t") # Keep time, pressure, temperature
  df$time <- matlab2posix(df$tim)
  return(df)
}

# Batch Loading ctd
#D02 - Mad01-18 ctd044-061
#D03 - Mad19-20 ctd062-063

ctd_files <- list.files(
  "CTDdataMSM126andM209/msm_126_1/ctd/",
  pattern = "\\.ctd$",
  full.names = TRUE
)

ctd_files_sel <- ctd_files[
  as.numeric(sub(".*_(\\d{3})_1sec\\.ctd$", "\\1", ctd_files)) %in% 44:63
]

ctd_list <- lapply(ctd_files_sel, read_ctd_tim_p)

names(ctd_list) <- paste0("Mad", sprintf("%02d", seq_along(ctd_list)))


#Cleaning steps
#get the names of all files before applying both cleaning functions to all files

ctd_list_clean <- lapply(ctd_list, function(df) {
  df %>%
    clean_start_depth("p") %>%
    clean_end_ascent("p")
})

# Bind all casts
ctd_all <- bind_rows(ctd_list_clean, .id = "dive") %>%
  mutate(
    time = lubridate::as_datetime(time),
    date = as.Date(time)) %>%
  select(dive, time, p, ctd_temp = t)

#.---------------------------------------------------------------------------------------
# 2.4. Read CTD data for D04 ------------------------------------------------------------
#.---------------------------------------------------------------------------------------

# 1. Upload and processing of ctd data for the 3 dives of D4 ---------------------------------------------

#Mean temperature and depth for each timestamp because the ctd takes 3 measures every 0.5 s
#. AR1 ----------------------------------------------
ctd_raw1 <- read.csv("Scripts&Outputs/EnvDeep01_03_T7.3-01-03_vs_CTD_oom/clean_ctd/K2261300_processed_DESCENT_ASCENT/K2261300_processed_DESCENT_ASCENT.csv")

# Convert time
ctd_clean1 <- ctd_raw1 %>%
  mutate(
    time = ymd_hms(`DATETIME..YYYY.MM.DD.HH.MM.SS.` , tz = "UTC"),
    Station = recode(Station, "AR2" = "AR1")
  ) %>%
  # Keep only good quality
  filter(
    TEMP_QC == 1,
    PRES_QC == 1,
    PRES..dbar. >= 1.5
  ) %>%
  # Keep only relevant columns
  select(time, Station, pres = `PRES..dbar.`, temp = `TEMP...C.`) %>%
  # Aggregate replicates (same timestamp)
  group_by(time, Station) %>%
  summarise(
    p = mean(pres, na.rm = TRUE),
    ctd_temp = mean(temp, na.rm = TRUE),
    .groups = "drop"
  )


#. SN ----------------------------------------------
ctd_raw2 <- read.csv("Scripts&Outputs/EnvDeep01_03_T7.3-01-03_vs_CTD_oom/clean_ctd/K2261300_processed_DESCENT_ASCENT/K2261128_Sofia_500m_processed_DESCENT_ASCENT.csv")

ctd_clean2 <- ctd_raw2 %>%
  mutate(
    time = ymd_hms(`DATETIME..YYYY.MM.DD.HH.MM.SS.`, tz = "UTC"),
    Station = as.character(Station),
    Station = recode(Station, "CL3" = "SN", "PF1" = "SN")
  ) %>%
  filter(TEMP_QC == 1, PRES_QC == 1, PRES..dbar. >= 1.5) %>%
  select(time, Station, pres = `PRES..dbar.`, temp = `TEMP...C.`) %>%
  group_by(time, Station) %>%
  summarise(
    p = mean(pres, na.rm = TRUE),
    ctd_temp = mean(temp, na.rm = TRUE),
    .groups = "drop"
  )

#. SG3 ----------------------------------------------

ctd_raw3 <- read.csv("Scripts&Outputs/EnvDeep01_03_T7.3-01-03_vs_CTD_oom/clean_ctd/K2261300_processed_DESCENT_ASCENT/K2261400_processed_DESCENT_ASCENT.csv")

# Convert time
ctd_clean3 <- ctd_raw3 %>%
  mutate(
    time = ymd_hms(`DATETIME..YYYY.MM.DD.HH.MM.SS.` , tz = "UTC")
  ) %>%
  filter(
    TEMP_QC == 1,
    PRES_QC == 1,
    PRES..dbar. >= 1.5
  ) %>%
  select(time, Station, pres = `PRES..dbar.`, temp = `TEMP...C.`) %>%
  group_by(time, Station) %>%
  summarise(
    p = mean(pres, na.rm = TRUE),
    ctd_temp = mean(temp, na.rm = TRUE),
    .groups = "drop"
  )

# Combine files from D04
ctd_oom_all <- bind_rows(ctd_clean1, ctd_clean2, ctd_clean3) %>%
  rename(
    dive = Station)

#COMBINE ALL CTD CASTS
ctd_all <- bind_rows(ctd_all, ctd_oom_all)


#. ==============================================================================
# 4. Define descent / ascent per CTD cast --------------------------
#. ==============================================================================

ctd_phased <- ctd_all %>%
  group_by(dive) %>%
  mutate(
    p_max = max(p, na.rm = TRUE),
    t_pmax = time[p == p_max][1],
    phase = if_else(time <= t_pmax, "Descent", "Ascent")
  ) %>%
  ungroup()

#. ==============================================================================
# 5. Combine CTD data for all dives  -----------------------------------------
#. ==============================================================================
ctd_phased <- bind_rows(ctd_phased, rov01)  %>%
  select(dive, time, p, ctd_temp, phase)

#. ==============================================================================
# 6. Match CTD data (pressure and temperature) to logger's timestamps ------------
#. ==============================================================================

#Define time windows of each cast based on ctd
ctd_windows <- ctd_phased %>%
  group_by(dive) %>%
  summarise(
    time_start = min(time),
    time_end   = max(time),
    .groups = "drop")

#Create unique ids for dive x loggers (bc the same loggers are used in multiple dives)
pairs <- ctd_phased %>%
  distinct(dive) %>%
  crossing(loggers_all %>% distinct(logger))

match_ctd_logger <- function(dive_id, logger_id) {
  ctd_df <- ctd_phased %>%
    filter(dive == dive_id)
  win <- ctd_windows %>%
    filter(dive == dive_id)
  log_df <- loggers_all %>%
    filter(logger == logger_id)
  # restrict logger to CTD window (correct scope now)
  log_df <- log_df %>%
    filter(time >= win$time_start,
           time <= win$time_end)
  # nearest CTD match
  ctd_dt <- as.data.table(ctd_df)
  log_dt <- as.data.table(log_df)
  setkey(ctd_dt, time)
  setkey(log_dt, time)
  merged <- ctd_dt[log_dt, on = "time", roll = "nearest"]
  merged$dive <- dive_id
  merged$logger <- logger_id
  as.data.frame(merged)
}

paired_ctd <- mapply(
  FUN = match_ctd_logger,
  dive_id = pairs$dive,
  logger_id = pairs$logger,
  SIMPLIFY = FALSE
) %>%
  bind_rows()

#Check if times are correct
paired_summary <-paired_ctd %>%
  group_by(logger, dive) %>%
  summarise(
    start = min(time),
    end   = max(time),
    n     = n(),
    .groups = "drop"
  ) %>%
  arrange(logger, dive)


filename01 <- paste0("paired_all_", Sys.Date(), ".csv")
pathfile01 <- "C:/Users/asnog/Documents/PhD-MareMadeira/TemperatureLogger/Scripts&Outputs"
write.csv(paired_ctd, file = paste0(pathfile01, "/", filename01), row.names = FALSE)

#. ==============================================================================
# 7. Validation & Metadata metrics  -----------------------------------------
#. ==============================================================================

# Calculating DELTA T + Renaming variables for clarity
paired_all <- paired_ctd %>%
  rename(depth_m = p) %>%
  mutate(
    unique_id = paste(logger, dive, sep = "_"),
    delta_T = log_temp - ctd_temp
    
  )

# Rough thermocline estimation (per dive + logger + phase) ------------------------
thermocline_est <- ctd_phased %>%
  filter(p <= 200) %>%
  group_by(dive, phase) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    # rolling smoothing (your idea)
    temp_smooth = rollmean(ctd_temp, k = 7, fill = NA, align = "center"),
    p_smooth    = rollmean(p,        k = 7, fill = NA, align = "center")
  ) %>%
  mutate(
    # gradient on smoothed signal
    dT = c(NA, diff(temp_smooth)),
    dp = c(NA, diff(p_smooth)),
    dT_dz = ifelse(dp == 0, NA, dT / dp)
  ) %>%
  summarise(
    idx = which.max(abs(dT_dz)),
    thermocline_depth_m = p[idx],
    thermocline_time    = time[idx],
    max_gradient        = max(abs(dT_dz), na.rm = TRUE),
    .groups = "drop"
  )

# Phase-based performance metrics --------------------------------------
phase_summary <- paired_all %>%
  group_by(dive, phase, logger) %>%
  arrange(time, .by_group = TRUE) %>%
  summarise(
    # Time coverage of logger data
    time_start = min(time, na.rm = TRUE),
    time_end   = max(time, na.rm = TRUE),
    duration_s = as.numeric(difftime(time_end, time_start, units = "secs")),
    # Depth coverage of logger data
    depth_min_m = depth_m[which.min(time)],
    depth_max_m = depth_m[which.max(time)],
    depth_range_m = depth_max_m - depth_min_m,
    #Vertical Velocity
    vertical_speed_ms_1 = depth_range_m /duration_s,
    # Performance metrics
    bias = mean(delta_T, na.rm = TRUE),
    rmse = sqrt(mean(delta_T^2, na.rm = TRUE)),
    mae_phase  = mean(abs(delta_T), na.rm = TRUE),
    sd   = sd(delta_T, na.rm = TRUE),
    n    = n(),
    # Extremes
    max_abs_delta_T = max(abs(delta_T), na.rm = TRUE),
    time_at_max_abs_dT = time[which.max(abs(delta_T))],
    depth_at_max_abs_dT = depth_m[which.max(abs(delta_T))],
    .groups = "drop"
  )

phase_summary <- phase_summary %>%
  left_join(
    thermocline_est,
    by = c("dive", "phase")
  )

phase_summary <- phase_summary %>%
  group_by(dive, phase, logger) %>%
  mutate(
    lag_thermocline_s = as.numeric(
      difftime(time_at_max_abs_dT, thermocline_time, units = "secs") # timelag: logger vs thermocline event
    ),
    depth_lag_m = depth_at_max_abs_dT - thermocline_depth_m # depth lag
  ) 

#Could have been done sooner..
paired_all <- paired_all %>%
  mutate(
    logger = recode(logger,
                    "EnvLogger01" = "EnvTide01",
                    "EnvLogger02" = "EnvTide02",
                    "EnvLogger03" = "EnvTide03"
    )
  )

#Export
filename02 <- paste0( format(Sys.Date(), "%d-%m-%Y"), "_phase-summary-validation-metrics_all.csv")
pathfile02 <- "C:/Users/asnog/Documents/PhD-MareMadeira/TemperatureLogger/Scripts&Outputs/"
write.csv(phase_summary, file = paste0(pathfile02, "/", filename02), row.names = FALSE)


#. ==============================================================================
# 8. Overall logger performance plots  ------------------------------------------
#. ==============================================================================

#8.1. CORRELATION PLOT ----------------------------------------------

stats_by_logger_phase <- paired_all %>%
  group_by(logger, phase) %>%
  summarise(
    r2 = cor(ctd_temp, log_temp, use = "complete.obs")^2,
    bias = mean(delta_T, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "R²=", round(r2, 2), "\n",
      "\u03B4=", round(bias, 2), "°C\n",
      "n=", n
    )
  )

stats_by_logger_phase <- stats_by_logger_phase %>%
  mutate(
    x_pos = case_when(
      phase == "Descent" ~ 0.65,
      phase == "Ascent"  ~ 8.2,
      phase == "Bottom"  ~ 0.65
    ),
    y_pos = case_when(
      phase %in% c("Descent", "Ascent") ~ 19.75,
      phase == "Bottom" ~ 14
    )
  )

# ---- plot ----
correlation_plots <- ggplot(paired_all, aes(x = ctd_temp, y = log_temp, color = phase)) +
  geom_point(alpha = 0.2, size = 0.7) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "longdash", color = "#333333", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.85, alpha = 0.7) +
  scale_color_manual(values = c(
    "Descent" = "#8D1C06FF",
    "Ascent"  = "#2C3778FF",
    "Bottom"  = "#DABD61FF" #FDDFA4FF
  )) +
  facet_wrap(~logger) +
  coord_equal(xlim = c(0, 22.5), ylim = c(0, 22.5), expand = FALSE) +
  labs(
    title = "Logger Correlation: Logger vs. CTD",
    subtitle = "Dashed line represents 1:1 perfect agreement",
    x = "CTD Temperature (°C)",
    y = "Logger Temperature (°C)",
    color = "Dive Phase"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.margin = margin(t = 5, r = 5, b = 2, l = 1),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    panel.spacing.x = unit(1, "lines"),
    legend.key.size = unit(1.2, "cm"),
    legend.text = element_text(size = 16),
    legend.margin = margin(t = -6),
    axis.title.x.bottom = element_text(margin = margin(t = 10)),
    axis.title.y.left = element_text(margin = margin(r = 15))
    ) +
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
  geom_text(
    data = stats_by_logger_phase,
    aes(x = x_pos, y = y_pos, label = label, color = phase),
    hjust = 0, size = 4, fontface = "bold", lineheight = 0.9, inherit.aes = FALSE, show.legend = FALSE
  )

print(correlation_plots)


#8.2. RESIDUALS HISTOGRAM ---------------------------------------------------
# Calculate stats for histogram labels
stats_hist <- paired_all %>%
  group_by(logger, phase) %>%
  summarise(
    bias = mean(delta_T, na.rm = TRUE),
    sd   = sd(delta_T, na.rm = TRUE),
    n    = n(),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("\u03B4 = ", round(bias, 2), " (\u00B1", round(sd, 2), ")")
  ) %>%
  mutate( #Adjust label based on the values specific for this this dataset
    x_pos = -4.2, 
    y_pos = case_when(
      phase == "Descent" ~ 0.165,
      phase == "Bottom"  ~ 0.135,
      phase == "Ascent"  ~ 0.15
    )
  )


residual_histograms <- ggplot(paired_all, aes(x = delta_T, fill = phase)) +
  # This formula forces the sum to be per facet (PANEL)
  geom_histogram(
    aes(y = after_stat(count) / tapply(after_stat(count), PANEL, sum)[PANEL]), 
    bins = 75, 
    alpha = 0.6, 
    position = "identity", 
    color = "white", 
    linewidth = 0.2
  ) +
  # Zero-error reference line
  geom_vline(xintercept = 0, linetype = "longdash", color = "#333333", linewidth = 0.8) +
  # Colors matching your correlation plot
  scale_fill_manual(values = c(
    "Descent" = "#8D1C06FF",
    "Ascent"  = "#2C3778FF",
    "Bottom"  = "#DABD61FF"
  )) +
  # Formatting Y axis as Percent
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::percent_format(), breaks = seq(0, 0.15, by = 0.05)) +
  coord_cartesian(xlim = c(-4, 4), ylim=c(0, 0.17)) + 
  facet_wrap(~logger) +
  labs(
    title = "Distribution of Temperature Differences",
    subtitle = expression(Delta*T == T[logger] - T[ctd]),
    x = "Temperature Difference (°C)",
    y = "Frequency (%)",
    fill = "Dive Phase"
  ) +
  theme_bw(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
  ) +
  # Add the stats labels
  geom_text(
    data = stats_hist,
    aes(x = x_pos, y = y_pos, label = label, color = phase),
    hjust = 0, size = 4.5, fontface = "bold", inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Descent" = "#8D1C06FF",
    "Ascent"  = "#2C3778FF",
    "Bottom"  = "#DABD61FF"
  ))

print(residual_histograms)


#. ================================================================================================
# 9. Exploratory analysis for temperature gradient and speed --------------------------------------
#. ================================================================================================

# 1. Add instantaneous gradients to your paired data
paired_all <- paired_all %>%
  group_by(unique_id, phase) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    # Time step in seconds
    dt_sec = as.numeric(difftime(time, lag(time), units = "secs")), #calculates how many secs bet "now" and previous timestamp
    # Change in CTD temperature and pressure
    dTemp = ctd_temp - lag(ctd_temp), #calculates change in true (ctd) water temperature since the last reading
    dPress = depth_m - lag(depth_m), #calculates how many meters the logger moved since the last reading
    # Instantaneous Rates of Change
    rate_temp_change = dTemp / dt_sec,  # °C per second (Plot 1)
    vertical_speed   = dPress / dt_sec  # meters per second (Plot 2)
  ) %>%
  ungroup()

#. -----------------------------------------------------------------------------
# Plot 1 Error (Temp Diff) vs. Rate of Temperature Change °C per second
#. -----------------------------------------------------------------------------

# Sanity check - filter out NA rows, if any, created by the lag() function
paired_analysis <- paired_all %>% filter(!is.na(rate_temp_change))

# Calculate R2 for separately for EACH logger AND phase to see how well the linear model fits
stats_phase_r2 <- paired_analysis %>%
  group_by(logger, phase) %>%
  summarise(
    r_sq = cor(rate_temp_change, delta_T, use = "complete.obs")^2,
    # The slope is a proxy for the 'Time Constant'
    slope = coef(lm(delta_T ~ rate_temp_change))[2], 
    .groups = "drop"
  ) %>%
  mutate(label = paste0("R² = ", round(r_sq, 2)))%>%
  mutate( #Adjust label based on the values specific for this this dataset
    x_pos = case_when(
      phase == "Descent" ~ -0.2,
      phase == "Bottom"  ~ 0.05,
      phase == "Ascent"  ~ 0.07
    ), 
    y_pos = case_when(
      phase == "Descent" ~ 8,
      phase == "Bottom"  ~ 2.5,
      phase == "Ascent"  ~ -8
    )
  )

# 2. Plot with separate lines for Ascent and Descent
error_vs_gradient_phased <- ggplot(paired_analysis, 
                                   aes(x = rate_temp_change, y = delta_T, color = phase)) +
  # Small points in background
  geom_point(alpha = 0.15, size = 0.5) + 
  
  # SEPARATE regression lines for each phase
  geom_smooth(method = "lm", linewidth = 0.6, alpha = 0.2, se = FALSE) +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Add the R2 label to each facet
  geom_text(data = stats_phase_r2, aes(x = x_pos, y = y_pos, label = label), 
            hjust = 0, fontface = "bold", size = 4.5, show.legend = FALSE) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0.05,0.05)) +
  coord_cartesian(xlim = c(-0.2, 0.2), ylim=c(-10, 10)) +
  facet_wrap(~logger) +
  scale_color_manual(values = c("Descent" = "#8D1C06FF", "Ascent" = "#2C3778FF", "Bottom" = "#DABD61FF")) +
  labs(
    title = c("Influence of the Rate of Temperature Change (°C s-1) on Temperature Differences"),
    subtitle = "Separate regressions per logger and per phase",
    x = "Rate of Temp Change (°C/sec)",
    y = expression(Delta*T ~ (Logger - CTD))
  ) +
  theme_bw(base_size = 14) +
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
  theme(strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
  )

print(error_vs_gradient_phased)


#. -----------------------------------------------------------------------------
# Plot 2: Error (Temp Diff) vs. Vertical Speed (m/s)
#. -----------------------------------------------------------------------------

# Sanity check - filter out NA rows
paired_analysis <- paired_all %>% filter(!is.na(vertical_speed))

# Calculate R2 for separately for EACH logger AND phase
stats_speed_r2 <- paired_analysis %>%
  group_by(logger, phase) %>%
  summarise(
    r_sq = cor(vertical_speed, delta_T, use = "complete.obs")^2,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("R² = ", round(r_sq, 2))) %>%
  mutate( 
    # Adjusted X positions based on typical vertical speeds (m/s)
    # Descent is usually positive speed, Ascent is negative speed
    x_pos = case_when(
      phase == "Descent" ~ 0.8,
      phase == "Bottom"  ~ -1.3,
      phase == "Ascent"  ~ -1.8
    ), 
    y_pos = case_when(
      phase == "Descent" ~ 5.1,
      phase == "Bottom"  ~ 2.8,
      phase == "Ascent"  ~ -5.3
    )
  )

# Plot with separate lines for Ascent and Descent
error_vs_speed_phased <- ggplot(paired_analysis, 
                                aes(x = vertical_speed, y = delta_T, color = phase)) +
  # Small points in background
  geom_point(alpha = 0.15, size = 0.5) + 
  
  # SEPARATE regression lines for each phase
  geom_smooth(method = "lm", linewidth = 0.6, alpha = 0.2, se = FALSE) +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  
  # Add the R2 label to each facet
  geom_text(data = stats_speed_r2, aes(x = x_pos, y = y_pos, label = label), 
            hjust = 0, fontface = "bold", size = 4.5, show.legend = FALSE) +
  
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0.1,0.1)) +
  
  # Adjusted xlim for vertical speed (m/s) - tweak if your max speed is different!
  coord_cartesian(xlim = c(-1.5, 1.5), ylim=c(-6, 6)) +
  
  facet_wrap(~logger) +
  scale_color_manual(values = c("Descent" = "#8D1C06FF", "Ascent" = "#2C3778FF", "Bottom" = "#DABD61FF")) +
  labs(
    title = "Influence of Vertical Speed (m/s) on Temperature Differences",
    subtitle = "Separate regressions per logger and per phase",
    x = "Vertical Speed (m/s)",
    y = expression(Delta*T ~ (Logger - CTD))
  ) +
  theme_bw(base_size = 14) +
  guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
  theme(strip.background = element_rect(fill = "grey90"),
        legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
  )

print(error_vs_speed_phased)

#. ================================================================================================
# 10. Instrument Comparison - How comparable are the loggers when deployed together? -------------------------------------------
#. ================================================================================================

#1. Introduce deployment column for better filtering
paired_analysis <- paired_analysis %>%
  mutate(deployment = case_when(
    dive == "ROV Dive" ~ "Depl_1",
    dive %in% paste0("Mad", sprintf("%02d", 1:18)) ~ "Depl_2",
    dive %in% c("Mad19", "Mad20") ~ "Depl_3",
    dive %in% c("AR1", "SG3", "SN") ~ "Depl_4",
    TRUE ~ "Other"
  ))


#2.1 Side-by-Side Comparison for Deployment #3 (Deep01 vs Deep02 at 4000m)

# Filter for Deployment3 
depl3_data <- paired_analysis %>% filter(deployment == "Depl_3")

# Stats for Deployment 3
depl3_stats <- depl3_data %>%
  group_by(logger) %>%
  summarise(
    n = n(),
    mae = mean(abs(delta_T), na.rm = TRUE),
    sd = sd(abs(delta_T), na.rm = TRUE),
    max_error = max(abs(delta_T), na.rm = TRUE),
    .groups = "drop"
  )

print(depl3_stats)

# Wilcoxon rank sum test for Deployment 3 only - paired = FALSE (default)
wilcox_depl3 <- wilcox.test(abs(delta_T) ~ logger, data = depl3_data)
print(wilcox_depl3)



#2.2. Comparison for Deployment #4 ( 2 Deep01 vs 3 Tide at 200-500m)

# 1. Add model_type to main dataframe
paired_analysis <- paired_analysis %>%
  mutate(model_type = if_else(str_detect(logger, "EnvTide"), "EnvTide", "EnvDeep"))

# Filter for Deployment 4 
depl4_data <- paired_analysis %>% filter(deployment == "Depl_4")

# Statistics for Deployment 4 only
depl4_individual_stats <- depl4_data %>%
  group_by(logger) %>%
  summarise(
    n = n(),
    bias = mean(delta_T, na.rm = TRUE),
    mae = mean(abs(delta_T), na.rm = TRUE),
    max_error = max(abs(delta_T), na.rm = TRUE),
    sd = sd(abs(delta_T), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

# Stats by Logger Type for Deployment 4 only
depl4_type_stats <- depl4_data %>%
  group_by(model_type) %>%
  summarise(
    n = n(),
    bias = mean(delta_T, na.rm = TRUE),
    mae = mean(abs(delta_T), na.rm = TRUE),
    sd = sd(abs(delta_T), na.rm = TRUE),
    max_error = max(abs(delta_T), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(logger = model_type) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

depl4_comparison <- bind_rows(depl4_individual_stats, depl4_type_stats)

#2.2.1 TESTING LOGGER TYPE - Wilcoxon Rank Sum test
wilcox_depl4 <- wilcox.test(abs(delta_T) ~ model_type, 
                            data = depl4_data) #%>% mutate(model_type = if_else(str_detect(logger, "EnvTide"), "EnvTide", "EnvDeep")))
print(wilcox_depl4)

# Optional: Visualize the comparison to support the test result
ggplot(depl4_data, aes(x = model_type, y = abs(delta_T), fill = model_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 1.5)) +
  labs(
    title = "Error Distribution by Logger Type (Deployment 4)",
    subtitle = paste("Wilcoxon p-value:", format.pval(wilcox_depl4$p.value)),
    x = "Model Type",
    y = "|ΔT| (Absolute Error °C)"
  ) +
  theme_minimal()


#2.2.2 COMPARING EACH LOGGER INDIVIDUAL PERFORMANCE

#KRUSKAL-WALLIS (basically an extension of Wilcoxon Rank Sum test to more than 2 groups)
#test statistical difference between the medians of three or more independent groupS

#DUNN'S TEST - post-hoc test if Kruskal-Wallis detects differences
#performs pairwise comparisons between each independent groups to determine which are different.

# Kruskal-Wallis for units in Deployment 4 only
kw_depl4 <- kruskal.test(abs(delta_T) ~ logger, data = depl4_data)
print(kw_depl4)

# Post-hoc to see which specific units differ in this "fair fight"

dunn_depl4 <- dunn.test(abs(depl4_data$delta_T), 
                        depl4_data$logger, 
                        method="bh")
print(dunn_depl4)
# -------------------------------------------------------
#FYI overall mae values per logger
mae_summary <- paired_analysis %>%
  group_by(logger) %>%
  summarise(
    mae = mean(abs(delta_T), na.rm = TRUE),
    bias = mean(delta_T, na.rm = TRUE),
    rmse = sqrt(mean(delta_T^2, na.rm = TRUE)),
    sd   = sd(delta_T, na.rm = TRUE),
    n    = n()
  )

print(mae_summary)

mae_summary_model <- paired_analysis %>%
  group_by(model_type) %>%
  summarise(
    mae = mean(abs(delta_T), na.rm = TRUE),
    bias = mean(delta_T, na.rm = TRUE),
    rmse = sqrt(mean(delta_T^2, na.rm = TRUE)),
    sd   = sd(delta_T, na.rm = TRUE),
    n    = n()
  )

print(mae_summary_model)

#. ================================================================================================
# 10. Testing how deployment affects logger performance -------------------------------------------
#. ================================================================================================

# Using only EnvDeep01 since it is the only one used in all four configurations
# Does it perform better depending on conditions (ROV vs CTD; 200 m vs 4000 m)

deep01_comparison <- paired_analysis %>%
  filter(logger == "EnvDeep01") %>%
  group_by(deployment) %>%
  summarise(
    n = n(),
    mae = mean(abs(delta_T), na.rm = TRUE),
    max_error = max(abs(delta_T), na.rm = TRUE),
    sd   = sd(delta_T, na.rm = TRUE),
    .groups = "drop"
  )

print("--- EnvDeep01: Performance by Deployment Type ---")
print(deep01_comparison)

# Filter for EnvDeep01 across all deployments
deep01_all_depls <- paired_analysis %>% filter(logger == "EnvDeep01")

# Kruskal-Wallis: Did EnvDeep01's error vary significantly by deployment type?
kw_deep01 <- kruskal.test(abs(delta_T) ~ deployment, data = deep01_all_depls)
print(kw_deep01)

# Post-hoc Dunn Test to see WHICH deployments were different
dunn_deep01 <- dunn.test(abs(deep01_all_depls$delta_T), 
                         deep01_all_depls$deployment, 
                         method="bh")

#. ====================================================================================================================
# 11. BIN-GROUPING-BASED Vertical Plots For Multiple-Cast Deployments (D2-D4)------------------------------------------
#. ====================================================================================================================

make_bins <- function(df, depth_var = "depth_m") {
  df %>%
    mutate(
      p_bin = floor(.data[[depth_var]] / 10) * 10 + 5
    )
}

#. ---------------------------------------------------------------------------
# Depth vs Temperature Plot --------------------------------------------------
#. ---------------------------------------------------------------------------
plot_temp_profile <- function(depl_id) {
  data_sub <- paired_analysis %>% 
    filter(deployment == depl_id)
  # --- CTD ---
  ctd_data <- data_sub %>%
    make_bins("depth_m") %>%
    group_by(p_bin) %>%
    summarise(
      mean_val = mean(ctd_temp, na.rm = TRUE),
      sd_val   = sd(ctd_temp, na.rm = TRUE),
      category = "CTD",
      .groups = "drop"
    )
  # --- Logger ---
  logger_data <- data_sub %>%
    filter(phase %in% c("Descent", "Ascent")) %>%
    make_bins("depth_m") %>%
    group_by(p_bin, model_type, phase) %>%
    summarise(
      mean_val = mean(log_temp, na.rm = TRUE),
      sd_val   = sd(log_temp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      category = paste(model_type,
                       if_else(phase == "Descent", "Downcasts", "Upcasts")))
  plot_df <- bind_rows(ctd_data, logger_data) %>%
    group_by(category) %>%
    arrange(p_bin) %>%
    mutate(
      mean_val = zoo::rollmean(mean_val, k = 3, fill = NA, align = "center"),
      sd_val   = zoo::rollmean(sd_val,   k = 3, fill = NA, align = "center")
    ) %>%
    ungroup()
  base_colors <- c("CTD" = "#001D3E")
  logger_colors <- c(
    "EnvDeep Downcasts" = "#C6E309",
    "EnvDeep Upcasts"   = "#95bf25",
    "EnvTide Downcasts" = "#C57890",
    "EnvTide Upcasts"   = "#993357")
  ggplot(plot_df, aes(y = p_bin, x = mean_val, color = category, fill = category)) +
    geom_ribbon(aes(xmin = mean_val - sd_val, xmax = mean_val + sd_val),
                alpha = 0.2, color = NA) +
    geom_path(linewidth = 1.1) +
    scale_y_reverse(
      limits = c(max(data_sub$depth_m, na.rm = TRUE) + 10, 0),
      expand = c(0, 0)
    ) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_x_continuous(position = "top",expand = c(0, 0), limits = c(0, 22.5)) +
    scale_color_manual(values = c(base_colors, logger_colors)) +
    scale_fill_manual(values = c(base_colors, logger_colors)) +
    labs(title = paste("Temperature Profile –", depl_id), x = "Temperature (°C)", y = "Depth (m)") +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 28, margin = margin(b = 14)),
      axis.title.x.top = element_text(size = 28, margin = margin(b = 20)),
      legend.position = "bottom",
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 16),
      axis.title.y = element_text(size = 28, margin = margin(r = 20)),
      axis.text = element_text(size = 28),
      axis.line = element_line(color = "darkgrey", linewidth = 1),
      plot.margin = margin(5.5, 15.5, 5.5, 5.5) #t, r, b, l
    )
}

plot_temp_profile("Depl_2")
plot_temp_profile("Depl_3")
plot_temp_profile("Depl_4")


#. ----------------------------------------------------------------------------
#. Depth vs Delta-Temperature Plot --------------------------------------------
#. ----------------------------------------------------------------------------

plot_delta_profile <- function(depl_id) {
  data_sub <- paired_analysis %>% 
    filter(deployment == depl_id)
  # --- Delta ---
  delta_data <- data_sub %>%
    filter(phase %in% c("Descent", "Ascent")) %>%
    make_bins("depth_m") %>%
    group_by(p_bin, model_type, phase) %>%
    summarise(
      mean_val = mean(delta_T, na.rm = TRUE),
      sd_val   = sd(delta_T, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      category = paste(model_type,
                       if_else(phase == "Descent", "Downcasts", "Upcasts")))
  plot_df <- delta_data %>%
    group_by(category) %>%
    arrange(p_bin) %>%
    mutate(
      mean_val = zoo::rollmean(mean_val, k = 3, fill = NA, align = "center"),
      sd_val   = zoo::rollmean(sd_val,   k = 3, fill = NA, align = "center")
    ) %>%
    ungroup()
  logger_colors <- c(
    "EnvDeep Downcasts" = "#C6E309",
    "EnvDeep Upcasts"   = "#95bf25",
    "EnvTide Downcasts" = "#D18E95",
    "EnvTide Upcasts"   = "#772B4E" 
  )
  ggplot(plot_df, aes(y = p_bin, x = mean_val, color = category, fill = category)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "#333333") +
    geom_ribbon(aes(xmin = mean_val - sd_val, xmax = mean_val + sd_val), alpha = 0.4, color = NA) +
    geom_path(linewidth = 1.1) +
    scale_y_reverse(limits = c(max(data_sub$depth_m, na.rm = TRUE) + 10, 0), expand = c(0, 0)) +
    scale_x_continuous(position = "top", expand = c(0, 0), limits = c(-4.5, 4.5)) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_color_manual(values = logger_colors) +
    scale_fill_manual(values = logger_colors) +
    labs(title = paste("ΔT Profile –", depl_id), x = expression(Delta * "T (Logger - CTD, °C)"), y = "Depth (m)"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 28, margin = margin(b = 14)),
      axis.title.x.top = element_text(size = 28, margin = margin(b = 20)),
      legend.position = "bottom",
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 28, margin = margin(r = 20)),
      axis.text = element_text(size = 28),
      axis.line = element_line(color = "darkgrey", linewidth = 1),
      plot.margin = margin(5.5, 15.5, 5.5, 5.5)
    )
}

#Call plots
plot_delta_profile("Depl_2")
plot_delta_profile("Depl_3")
plot_delta_profile("Depl_4")


#. ====================================================================================================================
# 12. Vertical Plots For ROV Deployment (D1) --------------------------------------------------------------------------
#. ====================================================================================================================

rov_data <- paired_analysis %>% filter(deployment == "Depl_1")

# Time-series plot -------------------------------------------------------------

plot_timeseries_rov <- function(df) {
  # Ranges
  temp_range  <- range(c(df$ctd_temp, df$log_temp), na.rm = TRUE)
  depth_range <- range(df$depth_m, na.rm = TRUE)
  # Scale + invert depth
  df <- df %>%
    mutate(depth_scaled = temp_range[1] +
        ((max(depth_m) - depth_m) / (depth_range[2] - depth_range[1])) *
        (temp_range[2] - temp_range[1]))
  ggplot(df, aes(x = time)) +
    # Depth
    geom_line(aes(y = depth_scaled), colour = "#6E6D70", linewidth = 1.2) +
    # CTD
    geom_line(aes(y = ctd_temp, colour = "CTD"), linewidth = 1.3, alpha = 0.9) +
    # Logger
    geom_point(aes(y = log_temp, colour = phase), size = 1.7, alpha = 0.6) +
    scale_colour_manual(
      values = c(
          "CTD"     = "#01244D",
          "Descent" = "#419F81",
          "Ascent"  = "#9B9EF5", ##542344
          "Bottom"  = "#dfa65b"),
      breaks = c("CTD", "Descent", "Bottom", "Ascent")) +
    scale_y_continuous(name = "Temperature (°C)",
      sec.axis = sec_axis(
        ~ depth_range[2] - (
          (. - temp_range[1]) / (temp_range[2] - temp_range[1]) *
            (depth_range[2] - depth_range[1])),
        name = "Depth (m)")) +
    guides(colour = guide_legend(override.aes = list(size = 6, alpha =1))) +
    labs(x = "Time (HH:MM)", title = "ROV Dive – Temperature & Depth Time Series") +
    theme_minimal(base_size = 20) +
    scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hours") +
    theme(
      axis.line = element_line(color = "#353535"),
      legend.position = "bottom",
      panel.grid.major.x = element_line(color = "grey80", linewidth = 0.5),
      axis.title.x = element_text(margin = margin(t = 15)),
      axis.title.y.left = element_text(margin = margin(r = 10)),
      axis.title.y.right = element_text(color = "#6E6D70", margin = margin(l = 10)),
      axis.text.y.right  = element_text(color = "#6E6D70")
    )
}

plot_timeseries_rov(rov_data)


#. ------------------------------------------------------------------------------------------
# Depth vs Temperature Plot -----------------------------------------------------------------
#. ------------------------------------------------------------------------------------------

plot_vertical_profile_rov <- function(df) {
  ggplot(df) +
    geom_point(aes(x = ctd_temp, y = depth_m, colour = "CTD"), size = 1, alpha=0.8) +
    geom_point(aes(x = log_temp, y = depth_m, colour = phase), alpha = 0.5, size = 1.5) +
    scale_y_reverse(limits = c(max(df$depth_m, na.rm = TRUE) + 15, 0), expand = c(0, 0)) +
    scale_x_continuous(position = "top", expand = c(0, 0), limits = c(0, 22.5)) +
    guides(colour = guide_legend(override.aes = list(size = 6, alpha =1, nrow = 4, byrow = TRUE))) +
    scale_colour_manual(
      values = c(
        "CTD"     = "#01244D",
        "Descent" = "#419F81",
        "Ascent"  = "#9B9EF5",
        "Bottom"  = "#dfa65b"
      ),  breaks = c("CTD", "Descent", "Bottom", "Ascent")) +
    labs(x = "Temperature (°C)",y = "Depth (m)",title = "ROV Vertical Temperature Profile") +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 28, margin = margin(b = 14)),
      axis.title.x.top = element_text(size = 28, margin = margin(b = 20)),
      legend.position = "right",
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 28, margin = margin(r = 25)),
      axis.title.x = element_text(size = 28, margin = margin(b = 15)),
      axis.text = element_text(size = 28),
      axis.line = element_line(color = "darkgrey", linewidth = 1),
      plot.margin = margin(5.5, 15.5, 5.5, 5.5))
}

plot_vertical_profile_rov (rov_data)


#. -------------------------------------------------------------------------------------------
# Depth vs ΔTemperature Plot (faceted by phase) ----------------------------------------------
#. -------------------------------------------------------------------------------------------

plot_delta_profile_rov <- function(df) {
  df <- df %>%
    mutate(
      phase = factor(phase, levels = c("Descent", "Bottom", "Ascent")))
  ggplot(df, aes(x = delta_T, y = depth_m, colour = phase)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_y_reverse(limits = c(max(df$depth_m, na.rm = TRUE) + 20, 0), expand = c(0, 0)) +
    scale_x_continuous(position = "bottom", limits = c(-4, 4), expand = c(0, 0), name = expression(Delta * "T (Logger - CTD, °C)", expand = c(0, 0))) +
    scale_colour_manual(
      values = c(
        "Descent" = "#419F81",
        "Ascent"  = "#9B9EF5",
        "Bottom"  = "#dfa65b"),
      breaks = c("Descent", "Bottom", "Ascent")) +
    facet_wrap(~phase, nrow = 1, strip.position = "top") +
    labs(y = "Depth (m)", title = "ΔT Profile by Phase (ROV)") +
    theme_minimal(base_size = 20) +
    theme(
      strip.placement = "outside",
      strip.text = element_text(size = 20, face ="bold", margin = margin(b = 15)),
      legend.position = "none",
      axis.title.y = element_text(size = 20, margin = margin(r = 20)),
      axis.title.x = element_text(size = 20, margin = margin(t = 120)),
      panel.spacing = unit(3, "lines"),
      axis.text = element_text(size = 20),
      axis.line = element_line(color = "#353535", linewidth = 0.8)
    )
}


plot_delta_profile_rov (rov_data)


#. ====================================================================================================================
# 13. Exploratory stats by LOGGER MODEL --------------------------------------------------------------------------
#. ====================================================================================================================

#Prepare data
#Make model a factor and order it:
paired_model <- paired_analysis %>%
  mutate(
    model_type = if_else(str_detect(logger, "EnvTide"), "EnvTide", "EnvDeep"),
    model_type = factor(model_type, levels = c("EnvDeep", "EnvTide"))
  )

phase_cols <- c(
  "Descent" = "#8D1C06FF",
  "Ascent"  = "#2C3778FF",
  "Bottom"  = "#DABD61FF"
)

# =============================================================================
# HISTOGRAM BY MODEL =========================================================
# =============================================================================

#Stats
stats_hist_model <- paired_model %>%
  group_by(model_type, phase) %>%
  summarise(
    bias = mean(delta_T, na.rm = TRUE),
    sd   = sd(delta_T, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("\u03B4 = ", round(bias, 2), " (\u00B1", round(sd, 2), ")"))

#plot
make_histogram <- function(model_name, xlim = c(-4, 4), ylim = c(0, 0.17)){
  x_span <- xlim[2] - xlim[1]
  y_span <- ylim[2] - ylim[1]
  text_data <- stats_hist_model %>%
    filter(model_type == model_name) %>%
    mutate(
      x_pos = xlim[1] - 1.6 * y_span,
      y_pos = case_when(
        phase == "Descent" ~ ylim[2] - 0.05 * y_span,
        phase == "Bottom"  ~ ylim[2] - 0.26 * y_span,
        phase == "Ascent"  ~ ylim[2] - 0.16 * y_span
      )
    )
  ggplot(paired_model %>% filter(model_type == model_name), aes(x = delta_T, fill = phase)) +
    geom_histogram(
      aes(y = after_stat(count) / sum(after_stat(count))),
      bins = 65, alpha = 0.6, position = "identity", color = "white", linewidth = 0.2
    ) +
    geom_vline(xintercept = 0, alpha = 0.6, linetype = "longdash", color = "#333333", linewidth = 0.6) +
    geom_text(
      data = text_data,
      aes(x = x_pos, y = y_pos, label = label, color = phase),
      inherit.aes = FALSE, hjust = 0, size = 4.5, fontface = "bold", show.legend = FALSE
    ) +
    scale_fill_manual(values = phase_cols) +
    scale_color_manual(values = phase_cols) +
    scale_y_continuous(labels = scales::percent_format(), breaks = seq(0, ylim[2], by = 0.05)) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(
      title = paste(model_name, "Logger - \u0394T Histogram"),
      x = expression(Delta*T~"(°C)"),
      y = "Frequency (%)"
    ) +
    guides(color = guide_legend(title = "Phase"), fill = guide_legend(title = "Phase")) +
    theme_bw(base_size = 16) +
    theme(strip.background = element_rect(fill = "grey90"))
}

# =============================================================================
# GRADIENT STATS
# =============================================================================
stats_grad_model <- paired_model %>%
  filter(!is.na(rate_temp_change)) %>%
  group_by(model_type, phase) %>%
  summarise(
    r_sq = cor(rate_temp_change, delta_T, use = "complete.obs")^2,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("R² = ", round(r_sq, 2)))

# =============================================================================
# GRADIENT FUNCTION
# =============================================================================
make_gradient <- function(model_name, xlim = c(-0.2, 0.2), ylim = c(-10, 10)){
  x_span <- xlim[2] - xlim[1]
  y_span <- ylim[2] - ylim[1]
  
  # Custom positional mappings for temperature change rates
  text_data <- stats_grad_model %>%
    filter(model_type == model_name) %>%
    mutate(
      x_pos = case_when(
        phase == "Bottom"  ~ xlim[1] + 0.75 * x_span,   # Center/Right
        phase == "Descent" ~ xlim[1] + 0.03 * x_span,   # Top-Left
        phase == "Ascent"  ~ xlim[2] - 0.02 * x_span    # Bottom-Right
      ),
      y_pos = case_when(
        phase == "Bottom"  ~ ylim[2] - 0.40 * y_span,   # Center/Right
        phase == "Descent" ~ ylim[2] - 0.02 * y_span,   # Top-Left
        phase == "Ascent"  ~ ylim[1] + 0.02 * y_span    # Bottom-Right
      ),
      hjust = case_when(
        phase == "Bottom"  ~ 0.5,
        phase == "Descent" ~ 0,
        phase == "Ascent"  ~ 1
      )
    )
  # Arrange data so "Bottom" rows are at the end (plotted last = on top)
  plot_data <- paired_model %>% 
    filter(model_type == model_name, !is.na(rate_temp_change)) %>%
    arrange(factor(phase, levels = c("Descent", "Ascent", "Bottom")))
  
  ggplot(
    paired_model %>% filter(model_type == model_name, !is.na(rate_temp_change)),
    aes(x = rate_temp_change, y = delta_T, color = phase)
  ) +
    geom_point(alpha = 0.25, size = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_text(
      data = text_data,
      aes(x = x_pos, y = y_pos, label = label, color = phase, hjust = hjust),
      inherit.aes = FALSE, size = 4.5, fontface = "bold", show.legend = FALSE
    ) +
    scale_color_manual(values = phase_cols) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(
      title = paste(model_name, "- Temp Change Rate"),
      x = expression(paste("Rate of Temp Change (", degree,"C/s)")),
      y = expression(Delta*T)
    ) +
    guides(color = guide_legend(title = "Phase"), fill = guide_legend(title = "Phase")) +
    theme_bw(base_size = 16)
}

# =============================================================================
# SPEED STATS
# =============================================================================
stats_speed_model <- paired_model %>%
  filter(!is.na(vertical_speed)) %>%
  group_by(model_type, phase) %>%
  summarise(
    r_sq = cor(vertical_speed, delta_T, use = "complete.obs")^2,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("R² = ", round(r_sq, 2)))

# =============================================================================
# SPEED FUNCTION
# =============================================================================
make_speed <- function(model_name, xlim = c(-1.5, 1.5), ylim = c(-6, 6)){
  x_span <- xlim[2] - xlim[1]
  y_span <- ylim[2] - ylim[1]
  # Custom positional mappings for vertical profiles
  text_data <- stats_speed_model %>%
    filter(model_type == model_name) %>%
    mutate(
      x_pos = case_when(
        phase == "Descent" ~ xlim[2] - 0.03 * x_span,   # Almost Top-Right
        phase == "Bottom"  ~ xlim[1] + 0.20 * x_span,   # Center-Top
        phase == "Ascent"  ~ xlim[1] + 0.005 * x_span    # Almost Bottom-Left
      ),
      y_pos = case_when(
        phase == "Descent" ~ ylim[2] - 0.02 * y_span,   # Almost Top-right
        phase == "Bottom"  ~ ylim[2] - 0.25 * y_span,   # Center-Top
        phase == "Ascent"  ~ ylim[1] + 0.01 * y_span    # Almost Bottom-Left
      ),
      hjust = case_when(
        phase == "Descent" ~ 1,
        phase == "Bottom"  ~ 0.5,
        phase == "Ascent"  ~ 0
      )
    )
  # Arrange data so "Bottom" rows are at the end (plotted last = on top)
  plot_data <- paired_model %>% 
    filter(model_type == model_name, !is.na(vertical_speed)) %>%
    arrange(factor(phase, levels = c("Descent", "Ascent", "Bottom")))
  ggplot(
    paired_model %>% filter(model_type == model_name, !is.na(vertical_speed)),
    aes(x = vertical_speed, y = delta_T, color = phase)
  ) +
    geom_point(alpha = 0.25, size = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_text(
      data = text_data,
      aes(x = x_pos, y = y_pos, label = label, color = phase, hjust = hjust),
      inherit.aes = FALSE, size = 4.5, fontface = "bold", show.legend = FALSE
    ) +
    scale_color_manual(values = phase_cols) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(
      title = paste(model_name, " - Vertical Speed"),
      x = "Vertical Speed (m/s)",
      y = expression(Delta*T)
    ) +
    guides(color = guide_legend(title = "Phase"), fill = guide_legend(title = "Phase")) +
    theme_bw(base_size = 16)
}

# =============================================================================
# CREATE THE 6 PANELS
# =============================================================================
# ---- EnvDeep Panels ----
p_hist_deep  <- make_histogram("EnvDeep", xlim = c(-4.5, 4.5), ylim = c(0, 0.10))
p_grad_deep  <- make_gradient("EnvDeep",  xlim = c(-0.2, 0.2),   ylim = c(-10, 10)) 
p_speed_deep <- make_speed("EnvDeep",     xlim = c(-2, 2),       ylim = c(-5.5, 5.5))

# ---- EnvTide Panels ----
p_hist_Tide  <- make_histogram("EnvTide", xlim = c(-2, 2),       ylim = c(0, 0.10))
p_grad_Tide  <- make_gradient("EnvTide",  xlim = c(-0.07, 0.07), ylim = c(-2, 2))
p_speed_Tide <- make_speed("EnvTide",     xlim = c(-1.8, 1.8),   ylim = c(-2, 2))

# =============================================================================
# PATCHWORK LAYOUT COMPILATION
# =============================================================================
final_model_panel <- (
  p_hist_deep + p_grad_deep + p_speed_deep
) / (
  p_hist_Tide + p_grad_Tide + p_speed_Tide
) +
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "center",
    axis.title.x = element_text(margin = margin(t = 5)), 
    axis.title.y = element_text(margin = margin(r = 5))
  )

final_model_panel

