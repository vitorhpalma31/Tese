library(gridExtra)


###############################################################################
## 1. Set number of simulations
###############################################################################
nsim <- 1000  # or however many you want

###############################################################################
## 2. Prepare storage for results
###############################################################################
frozen_results <- vector("list", length = length(mcmc_round_results))

###############################################################################
## 3. Main loop over each MCMC in mcmc_round_results
###############################################################################
for (k in seq_along(mcmc_round_results)) {
  
  # Which round was this MCMC frozen at?
  frozen_round <- (total_rounds / 2) + k - 1
  
  cat("\n=== Using MCMC #", k, " (frozen at round ", frozen_round, ") ===\n")
  
  # Extract the MCMC object
  round_mcmc <- mcmc_round_results[[k]]
  
  # Convert MCMC object to matrix of posterior samples
  all_samples_frozen <- as.matrix(round_mcmc$samples)
  
  # Future rounds to simulate
  start_round <- frozen_round + 1
  end_round   <- total_rounds
  
  # We'll store simulation-wise results
  correct_percentages   <- numeric(nsim)
  all_sim_standings     <- vector("list", nsim)
  all_sim_predictions   <- vector("list", nsim)
  
  #############################################################################
  ## Build partial standings from real data up to 'frozen_round'
  #############################################################################
  
  cat("Building actual standings from real data up to round", frozen_round, "\n")
  
  actual_standings <- data.frame(
    Team   = list_teams,
    Points = 0,
    Wins   = 0,
    Draws  = 0,
    Losses = 0
  )
  
  past_data <- subset(data, Round <= frozen_round)
  
  # Update actual_standings from real results
  for (i in seq_len(nrow(past_data))) {
    row_match <- past_data[i, ]
    h_team  <- row_match$HomeTeam
    a_team  <- row_match$AwayTeam
    h_goals <- row_match$ObsHGoals
    a_goals <- row_match$ObsAGoals
    
    if (is.na(h_goals) || is.na(a_goals)) {
      next  # skip if any NA in the real data
    }
    
    if (h_goals > a_goals) {
      actual_standings[actual_standings$Team == h_team, "Wins"]   <-
        actual_standings[actual_standings$Team == h_team, "Wins"] + 1
      actual_standings[actual_standings$Team == a_team, "Losses"] <-
        actual_standings[actual_standings$Team == a_team, "Losses"] + 1
    } else if (h_goals < a_goals) {
      actual_standings[actual_standings$Team == a_team, "Wins"]   <-
        actual_standings[actual_standings$Team == a_team, "Wins"] + 1
      actual_standings[actual_standings$Team == h_team, "Losses"] <-
        actual_standings[actual_standings$Team == h_team, "Losses"] + 1
    } else {
      # draw
      actual_standings[actual_standings$Team == h_team, "Draws"]   <-
        actual_standings[actual_standings$Team == h_team, "Draws"] + 1
      actual_standings[actual_standings$Team == a_team, "Draws"]   <-
        actual_standings[actual_standings$Team == a_team, "Draws"] + 1
    }
  }
  
  # Compute points
  actual_standings$Points <- actual_standings$Wins * 3 + actual_standings$Draws
  cat("Partial standings after round", frozen_round, ":\n")
  print(actual_standings)
  
  #############################################################################
  ## 4. Predict future rounds
  #############################################################################
  
  cat("Predicting future rounds from", start_round, "to", end_round, "\n")
  
  for (sim in seq_len(nsim)) {
    if (sim %% 50 == 0) {
      cat("   Simulation", sim, "of", nsim, "for MCMC #", k, "\n")
    }
    
    predictions <- data.frame()
    
    for (round_number in start_round:end_round) {
      
      round_data <- data[data$Round == round_number, 
                         c("HomeTeam","AwayTeam","HTeamIdx","ATeamIdx",
                           "ObsHGoals","ObsAGoals")]
      if (nrow(round_data) == 0) {
        cat("No matches found for round", round_number, "\n")
        next
      }
      
      # We'll apply an "apply" approach for each row
      preds <- apply(round_data, 1, function(row) {
        
        h_team <- as.integer(row["HTeamIdx"])
        a_team <- as.integer(row["ATeamIdx"])
        
        num_simulations <- 1000
        num_goal_sim    <- 1000
        
        # Random sample from the MCMC
        samp_idx <- sample(seq_len(nrow(all_samples_frozen)), 
                           size = num_simulations, replace = TRUE)
        
        sampled_home  <- all_samples_frozen[samp_idx, "home"]
        sampled_att_h <- all_samples_frozen[samp_idx, sprintf("att[%d]", h_team)]
        sampled_def_a <- all_samples_frozen[samp_idx, sprintf("def[%d]", a_team)]
        sampled_att_a <- all_samples_frozen[samp_idx, sprintf("att[%d]", a_team)]
        sampled_def_h <- all_samples_frozen[samp_idx, sprintf("def[%d]", h_team)]
        sampled_p1    <- all_samples_frozen[samp_idx, sprintf("p1[%d]", h_team)]
        sampled_p2    <- all_samples_frozen[samp_idx, sprintf("p2[%d]", a_team)]
        
        mu1 <- exp(sampled_home + sampled_att_h + sampled_def_a)
        mu2 <- exp(sampled_att_a + sampled_def_h)
        
        h_goals_sim <- numeric(num_goal_sim)
        a_goals_sim <- numeric(num_goal_sim)
        
        for (i in seq_len(num_goal_sim)) {
          h_goals_sim[i] <- rZMPC(1, mu = mu1[i], p = sampled_p1[i])
          a_goals_sim[i] <- rZMPC(1, mu = mu2[i], p = sampled_p2[i])
        }
        
        # Use na.rm=TRUE to ignore invalid draws
        HWinVal  <- mean(h_goals_sim > a_goals_sim, na.rm = TRUE) 
        DrawVal  <- mean(h_goals_sim == a_goals_sim, na.rm = TRUE)
        AWinVal  <- mean(h_goals_sim < a_goals_sim, na.rm = TRUE)
        
        c(
          HWin       = HWinVal,
          Draw       = DrawVal,
          AWin       = AWinVal,
          MeanHGoals = mean(h_goals_sim, na.rm = TRUE),
          MeanAGoals = mean(a_goals_sim, na.rm = TRUE)
        )
      })
      
      # Threshold logic for predicted outcomes
      predicted_result <- ifelse(
        abs(preds["HWin", ] - preds["AWin", ]) < 0.03, 
        "Draw",
        ifelse(
          (preds["HWin", ] > preds["AWin", ] & 
             preds["HWin", ] - preds["Draw", ] >= 0.02), 
          "Home Win",
          ifelse(
            (preds["AWin", ] > preds["HWin", ] & 
               preds["AWin", ] - preds["Draw", ] >= 0.02), 
            "Away Win", 
            "Draw"
          )
        )
      )
      
      # Actual result
      actual_result <- ifelse(
        round_data$ObsHGoals > round_data$ObsAGoals, "Home Win",
        ifelse(round_data$ObsHGoals == round_data$ObsAGoals, "Draw", "Away Win")
      )
      
      round_results <- data.frame(
        Round = round_number,
        round_data[, c("HomeTeam","AwayTeam")],
        t(preds),
        PredictedResult = predicted_result,
        ActualResult    = actual_result
      )
      
      predictions <- rbind(predictions, round_results)
    }
    
    # Evaluate correctness
    correct_comp <- predictions$PredictedResult == predictions$ActualResult
    correct_percentages[sim] <- mean(correct_comp, na.rm=TRUE) * 100
    
    # -------------------------------------------------------------------------
    # Combine "actual_standings" with these predictions
    # -------------------------------------------------------------------------
    
    sim_standings <- actual_standings
    
    for (team in list_teams) {
      home_matches <- predictions[predictions$HomeTeam == team, ]
      away_matches <- predictions[predictions$AwayTeam == team, ]
      
      w <- sum(home_matches$PredictedResult == "Home Win", na.rm=TRUE) +
        sum(away_matches$PredictedResult == "Away Win", na.rm=TRUE)
      d <- sum(home_matches$PredictedResult == "Draw", na.rm=TRUE) +
        sum(away_matches$PredictedResult == "Draw", na.rm=TRUE)
      l <- sum(home_matches$PredictedResult == "Away Win", na.rm=TRUE) +
        sum(away_matches$PredictedResult == "Home Win", na.rm=TRUE)
      
      sim_standings[sim_standings$Team == team, "Wins"]   <-
        sim_standings[sim_standings$Team == team, "Wins"] + w
      sim_standings[sim_standings$Team == team, "Draws"]  <-
        sim_standings[sim_standings$Team == team, "Draws"] + d
      sim_standings[sim_standings$Team == team, "Losses"] <-
        sim_standings[sim_standings$Team == team, "Losses"] + l
    }
    
    sim_standings$Points <- sim_standings$Wins * 3 + sim_standings$Draws
    sim_standings <- sim_standings[order(-sim_standings$Points, sim_standings$Team), ]
    sim_standings$Ranking <- seq_len(nrow(sim_standings))
    sim_standings <- sim_standings[, c("Ranking","Team","Points","Wins","Draws","Losses")]
    
    all_sim_predictions[[sim]] <- predictions
    all_sim_standings[[sim]]   <- sim_standings
  }
  
  #############################################################################
  ## 5. After all nsim, summarize
  #############################################################################
  
  avg_correct <- mean(correct_percentages, na.rm=TRUE)
  
  cat("=> For MCMC #", k,
      "(frozen at round", frozen_round, "):\n",
      "   average % correct =", round(avg_correct, 2), "%\n")
  
  frozen_results[[k]] <- list(
    MCMC_index      = k,
    round_frozen    = frozen_round,
    correct_percent = correct_percentages,
    avg_correct     = avg_correct,
    standings       = all_sim_standings,
    predictions     = all_sim_predictions
  )
  
} # end MCMC loop

###############################################################################
## 6. Compute overall mean across MCMCs if desired
###############################################################################
overall_mean_correct <- mean(sapply(frozen_results, function(x) x$avg_correct))
cat("Overall mean of avg_correct across all MCMCs:", overall_mean_correct, "\n")

###############################################################################
###############################################################################
###############################################################################


################################################################################
## 1. Create lists to hold the average standings and average points (one per MCMC).
################################################################################

n_mcmc <- length(frozen_results)   # e.g., 19 if you have 19 MCMCs
avg_standings_list <- vector("list", length = n_mcmc)
avg_points_list <- vector("list", length = n_mcmc)  # New list for average points

################################################################################
## 2. Loop over each MCMC in frozen_results
################################################################################

for (k in seq_along(frozen_results)) {
  
  # Extract the list of standings data frames for MCMC #k
  standings_list <- frozen_results[[k]]$standings
  
  # Number of simulations for this MCMC
  nsim <- length(standings_list)
  
  # Teams and number of teams
  teams   <- standings_list[[1]]$Team  # Extract team names from the first data frame
  n_teams <- length(teams)
  
  # Initialize a frequency matrix: rows = teams, columns = rank positions (1..n_teams)
  rankFreq <- matrix(
    0, 
    nrow = n_teams, 
    ncol = n_teams,
    dimnames = list(teams, as.character(seq_len(n_teams)))
  )
  
  # Initialize a points vector: rows = teams, to store total points per team
  totalPoints <- setNames(rep(0, n_teams), teams)
  
  ##############################################################################
  ## 2a. Fill the matrix with counts of how many times each Team got each Rank
  ##     and sum the points for each Team
  ##############################################################################
  for (sim_df in standings_list) {
    # Each sim_df is a data frame with columns like Team, Ranking, Points, etc.
    for (row_idx in seq_len(nrow(sim_df))) {
      this_team  <- sim_df$Team[row_idx]
      this_rank  <- sim_df$Ranking[row_idx]
      this_points <- sim_df$Points[row_idx]
      
      rankFreq[this_team, as.character(this_rank)] <-
        rankFreq[this_team, as.character(this_rank)] + 1
      
      # Accumulate points for this team
      totalPoints[this_team] <- totalPoints[this_team] + this_points
    }
  }
  
  # Convert counts to frequencies
  rankFreq <- rankFreq / nsim
  
  # Calculate the average points for each team
  avgPoints <- totalPoints / nsim
  
  # Round to 4 decimals
  rankFreq <- round(rankFreq, 4)
  avgPoints <- round(avgPoints, 4)
  
  ##############################################################################
  ## 2b. Prepare the average standings data frame
  ##############################################################################
  
  # 'rankFreq' is a matrix with rownames = team names, colnames = "1","2",...
  df_k <- as.data.frame(rankFreq)            # Wide data frame
  df_k$Team <- rownames(df_k)                # Put team names into a column
  
  # Reorder columns to put "Team" first, followed by rank columns
  df_k <- df_k[, c("Team", as.character(seq_len(n_teams)))]
  
  # Store the standings data frame in avg_standings_list
  avg_standings_list[[k]] <- df_k
  
  ##############################################################################
  ## 2c. Prepare the average points data frame
  ##############################################################################
  
  # Create a data frame for average points
  points_df_k <- data.frame(
    Team = names(avgPoints),
    AvgPoints = avgPoints
  )
  
  # Store the points data frame in avg_points_list
  avg_points_list[[k]] <- points_df_k
}

################################################################################
## 3. (Optional) Name each element of the lists so you know which MCMC it belongs to
################################################################################

names(avg_standings_list) <- paste0("avg_standings_MCMC_", seq_len(n_mcmc))
names(avg_points_list) <- paste0("avg_points_MCMC_", seq_len(n_mcmc))



################################################################################
################################################################################
################################################################################
# Step 1: Function to Calculate Observed Standings
calculate_observed_standings <- function(data_subset) {
  # Initialize the standings table
  standings <- data.frame(
    Team = unique(c(data_subset$HomeTeam, data_subset$AwayTeam)),
    Points = 0, Wins = 0, Draws = 0, Losses = 0,
    stringsAsFactors = FALSE
  )
  
  # Iterate over each team to calculate standings
  for (team in standings$Team) {
    # Count wins, draws, and losses for home and away games
    home_wins <- sum(data_subset$HomeTeam == team & data_subset$Result == "HomeWin")
    away_wins <- sum(data_subset$AwayTeam == team & data_subset$Result == "AwayWin")
    draws_as_home <- sum(data_subset$HomeTeam == team & data_subset$Result == "Draw")
    draws_as_away <- sum(data_subset$AwayTeam == team & data_subset$Result == "Draw")
    total_games_played <- sum(data_subset$HomeTeam == team | data_subset$AwayTeam == team)
    
    # Update the standings
    standings[standings$Team == team, ] <- data.frame(
      Team = team,
      Points = (home_wins + away_wins) * 3 + (draws_as_home + draws_as_away),
      Wins = home_wins + away_wins,
      Draws = draws_as_home + draws_as_away,
      Losses = total_games_played - (home_wins + away_wins + draws_as_home + draws_as_away)
    )
  }
  
  # Order teams by points
  standings <- standings[order(-standings$Points), ]
  return(standings)
}

# Step 2: Observed Standings for Different Parts of the Season
# Full observed standings for the entire season
full_observed_standings <- calculate_observed_standings(data)

# Observed standings for the first half of the season
observed_data_1_half <- data %>% filter(Round <= nteams - 1)
observed_standings_1_half <- calculate_observed_standings(observed_data_1_half)

# Observed standings for the second half of the season
observed_data_2_half <- data %>% filter(Round > nteams - 1)
observed_standings_2_half <- calculate_observed_standings(observed_data_2_half)
################################################################################
## Done! You now have 19 data frames inside avg_standings_list.
## Each data frame has columns: "Team", "1", "2", ..., "n_teams"
## with the probabilities of finishing in each rank.
################################################################################
################################################################################
##############################################################################
## Example: Compute probabilities for Top 4, Top 6, Top 7, and Bottom 3
##          for each MCMC in 'avg_standings_list'.
##############################################################################

# 1) Create a parent list 'top' with sub-lists for each category we want
top <- list(
  Top4     = vector("list", length = length(avg_standings_list)),
  Top6     = vector("list", length = length(avg_standings_list)),
  Top7     = vector("list", length = length(avg_standings_list)),
  Bottom3  = vector("list", length = length(avg_standings_list))
)

# 2) Loop over each MCMC (each element of 'avg_standings_list')
for (k in seq_along(avg_standings_list)) {
  
  # Extract the k-th MCMC's average standings data frame
  df_k <- avg_standings_list[[k]]
  
  # df_k has columns: "Team", "1", "2", ..., "20" (for a 20-team league)
  
  ## A) Compute Top 4 (ranks 1..4)
  prob_top4 <- rowSums(df_k[, c("1","2","3","4")])
  top4_df <- data.frame(
    Team            = df_k$Team,
    ProbabilityTop4 = prob_top4
  )
  # Filter out teams with zero probability
  top4_df <- subset(top4_df, ProbabilityTop4 > 0)
  
  ## B) Compute Top 6 (ranks 1..6)
  prob_top6 <- rowSums(df_k[, c("1","2","3","4","5","6")])
  top6_df <- data.frame(
    Team            = df_k$Team,
    ProbabilityTop6 = prob_top6
  )
  top6_df <- subset(top6_df, ProbabilityTop6 > 0)
  
  ## C) Compute Top 7 (ranks 1..7)
  prob_top7 <- rowSums(df_k[, c("1","2","3","4","5","6","7")])
  top7_df <- data.frame(
    Team            = df_k$Team,
    ProbabilityTop7 = prob_top7
  )
  top7_df <- subset(top7_df, ProbabilityTop7 > 0)
  
  ## D) Compute Bottom 3 
  prob_bottom3 <- rowSums(df_k[, tail(colnames(df_k), 3)])
  bottom3_df <- data.frame(
    Team              = df_k$Team,
    ProbabilityBottom3 = prob_bottom3
  )
  bottom3_df <- subset(bottom3_df, ProbabilityBottom3 > 0)
  
  # 3) Store each data frame in the respective sub-list of 'top'
  top$Top4[[k]]     <- top4_df
  top$Top6[[k]]     <- top6_df
  top$Top7[[k]]     <- top7_df
  top$Bottom3[[k]]  <- bottom3_df
}

# 4) (Optional) label each sub-list element by MCMC index
names(top$Top4)     <- paste0("MCMC_", seq_along(avg_standings_list))
names(top$Top6)     <- paste0("MCMC_", seq_along(avg_standings_list))
names(top$Top7)     <- paste0("MCMC_", seq_along(avg_standings_list))
names(top$Bottom3)  <- paste0("MCMC_", seq_along(avg_standings_list))

##############################################################################
## Done: 'top' is now a list with sub-lists:
##   top$Top4[[k]], top$Top6[[k]], top$Top7[[k]], top$Bottom3[[k]]
## Each element is a data frame of teams with >0% chance in that category.
##############################################################################
# Criar um vetor com os valores de avg_correct extraídos de frozen_results
avg_correct_values <- sapply(frozen_results, function(x) x$avg_correct)
mean(avg_correct_values)
# Gerar o gráfico
plot(avg_correct_values,
     type = "o",                # Tipo de gráfico com linhas e pontos
     col = "blue",              # Cor da linha e dos pontos
     lwd = 2,                   # Espessura da linha
     pch = 16,                  # Tipo de ponto preenchido
     xlab = "MCMC",             # Rótulo do eixo X
     ylab = "Correct Predictions %", # Rótulo do eixo Y
     main = "Correct Predictions for each MCMC", # Título do gráfico
     cex.lab = 1.3,             # Tamanho dos rótulos dos eixos
     cex.main = 1.6,            # Tamanho do título
     cex.axis = 1.1             # Tamanho dos números nos eixos
)
# Adicionar uma grade ao gráfico
grid()



