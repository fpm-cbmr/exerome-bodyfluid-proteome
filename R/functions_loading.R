# Define the mixed.anova function with covariates
mixed.anova <- function(dat, expo, outc, formel, covariates = NULL) {
    ## 'dat'  -- data frame containing all needed values
    ## 'expo' -- vector of exposure of interest
    ## 'outc' -- vector of outcomes of interest
    ## 'formel' -- string of the right-hand side of the formula (confounder and random effects)
    ## 'covariates' -- vector of covariates to be included in the model (default is NULL)

    ## load packages for mixed models
    library("lmerTest")

    ## Create a string of covariates if they are provided
    covariate_str <- if (!is.null(covariates)) {
        paste(covariates, collapse = " + ")
    } else {
        ""
    }

    ## run across all exposures and all outcomes of interest
    res <- lapply(outc, function(i) {
        ## loop over all exposures of interest
        tmp <- lapply(expo, function(j) {
            ## create formula with backticks around outcome and exposure
            ff <- if (covariate_str != "") {
                paste0("`", i, "` ~ `", j, "` + ", covariate_str, formel)
            } else {
                paste0("`", i, "` ~ `", j, "`", formel)
            }
            ## run the model
            ml <- lmer(as.formula(ff), data = dat)
            ff <- summary(ml)$coefficients
            ## return needed information
            ff <- data.frame(outcome = i, exposure = j, ind = row.names(ff), beta = ff[,1], se = ff[,2], pval = ff[,5])
            ## reshape the data
            ff <- reshape(ff, idvar = c("outcome", "exposure"), timevar = "ind", direction = "wide")
            ## exchange names to allow combination afterwards
            names(ff) <- gsub(j, "exposure", names(ff), fixed = TRUE)
            ## add number of observations
            ff$n <- nrow(na.omit(dat[, c(expo, i)]))
            ## add ANOVA p-values
            ml <- anova(ml)
            ff[, paste0("pval.aov.", rownames(ml))] <- ml$`Pr(>F)`
            return(ff)
        })
        ## combine into one data frame
        tmp <- do.call(rbind, tmp)
        return(tmp)
    })
    ## combine into one data frame
    res <- do.call(rbind, res)
    ## get only information needed
    res <- res[, c("outcome", "n", grep("exposure|aov", names(res), value = TRUE))]
    return(res)
}

# Define a function to perform mixed ANOVA in parallel
mixed_anova_parallel <- function(dat, expo, outc, formel, covariates = NULL) {
    res <- foreach(i = seq_along(outc), .combine = rbind, .packages = c("lmerTest", "dplyr", "reshape2"), .export = c("mixed.anova")) %dopar% {
        mixed.anova(dat, expo, outc[i], formel, covariates)  # Using mixed.anova function
    }
    return(res)
}


# group anova with handling for special characters in variable names
mixed.anova.group <- function(dat, expo, outc, formel, group_var, covariates = NULL) {
    library("lmerTest")

    # Create a string of covariates if they are provided, with each covariate wrapped in backticks
    covariate_str <- if (!is.null(covariates)) {
        paste(paste0("`", covariates, "`"), collapse = " + ")
    } else {
        ""
    }

    # Run across all exposures and all outcomes of interest
    res <- lapply(outc, function(i) {
        # Wrap outcome in backticks
        outcome <- paste0("`", i, "`")

        tmp <- lapply(expo, function(j) {
            # Wrap exposure and group_var in backticks
            exposure <- paste0("`", j, "`")
            group <- paste0("`", group_var, "`")

            # Construct the formula string with backticks around variables
            ff <- if (covariate_str != "") {
                paste0(outcome, " ~ ", exposure, " * ", group, " + ", covariate_str, formel)
            } else {
                paste0(outcome, " ~ ", exposure, " * ", group, formel)
            }

            # Run the model
            ml <- lmer(as.formula(ff), data = dat)
            ff <- summary(ml)$coefficients

            # Return needed information
            ff <- data.frame(outcome = i, exposure = j, ind = row.names(ff),
                             beta = ff[,1], se = ff[,2], pval = ff[,5])

            # Reshape the data
            ff <- reshape(ff, idvar = c("outcome", "exposure"), timevar = "ind", direction = "wide")
            # Exchange names to allow combination afterwards
            names(ff) <- gsub(j, "exposure", names(ff), fixed = TRUE)

            # Add number of observations
            ff$n <- nrow(na.omit(dat[, c(expo, i)]))

            # Add ANOVA p-values
            ml <- anova(ml)
            ff[, paste0("pval.aov.", rownames(ml))] <- ml$`Pr(>F)`  # Wrap `Pr(>F)` in backticks

            return(ff)
        })

        # Combine into one data frame
        tmp <- do.call(rbind, tmp)
        return(tmp)
    })

    # Combine all results into one data frame
    res <- do.call(rbind, res)
    # Get only information needed
    res <- res[, c("outcome", "n", grep("exposure|aov", names(res), value = TRUE))]
    return(res)
}

# Define the parallel function
mixed_anova_parallel_group <- function(dat, expo, outc, formel, group_var, covariates = NULL) {
    # Parallel execution with only the function exported
    res <- foreach(i = seq_along(outc), .combine = rbind, .packages = c("lmerTest", "dplyr", "reshape2"),
                   .export = "mixed.anova.group") %dopar% {
                       mixed.anova.group(dat, expo, outc[i], formel, group_var, covariates)
                   }
    return(res)
}

# function to run linear model for group seperately
run_model_for_group <- function(data, group_name, outc, formel, covariates = NULL) {
    # Subset the data for the specified group
    data_group <- data %>% filter(group == group_name)

    # Ensure `time` is a factor with `-1` as the baseline within this subset
    data_group <- data_group %>% mutate(time = factor(time, levels = c(-1, 0, 0.5, 24)))

    # Run the mixed ANOVA function on this subset
    res_group <- mixed_anova_parallel(
        dat = data_group,
        expo = "time",
        outc = outc,
        formel = formel,
        covariates = covariates
    )

    # Add the group label to the results for identification
    res_group <- res_group %>% mutate(group = group_name)

    return(res_group)
}

# function to run linear model for group seperately
run_model_for_sex <- function(data, sex_label, outc, formel, covariates = NULL) {
    # Subset the data for the specified group
    data_sex <- data %>% filter(sex == sex_label)

    # Ensure `time` is a factor with `-1` as the baseline within this subset
    data_sex <- data_sex %>% mutate(time = factor(time, levels = c(-1, 0, 0.5, 24)))

    # Run the mixed ANOVA function on this subset
    res_sex <- mixed_anova_parallel(
        dat = data_sex,
        expo = "time",
        outc = outc,
        formel = formel,
        covariates = covariates
    )

    # Add the group label to the results for identification
    res_sex <- res_sex %>% mutate(sex = sex_label)

    return(res_sex)
}


##
#Define the three-way ANOVA function with covariates
# Parallelized function to run three-way ANOVA in parallel
# Group-specific three-way ANOVA function to handle group-level analysis
mixed.anova.three_way_group <- function(dat, expo, outc, formel, group_var, covariates = NULL) {
    library("lmerTest")

    # Format covariates string
    covariate_str <- if (!is.null(covariates)) {
        paste(paste0("`", covariates, "`"), collapse = " + ")
    } else {
        ""
    }

    # Run three-way ANOVA across outcomes within each group
    res <- lapply(outc, function(i) {
        tmp <- lapply(expo, function(j) {
            # Include three-way interaction in formula: (group * sex * time)
            ff <- if (covariate_str != "") {
                paste0("`", i, "` ~ `", j, "` * ", group_var, " * sex * time + ", covariate_str, formel)
            } else {
                paste0("`", i, "` ~ `", j, "` * ", group_var, " * sex * time", formel)
            }

            # Run the model
            ml <- lmer(as.formula(ff), data = dat)
            summary_ml <- summary(ml)$coefficients

            # Create results data frame
            result <- data.frame(
                outcome = i,
                exposure = j,
                ind = row.names(summary_ml),
                beta = summary_ml[, 1],
                se = summary_ml[, 2],
                pval = summary_ml[, 5]
            )

            # Reshape for easy combination and add ANOVA p-values
            result <- reshape(result, idvar = c("outcome", "exposure"), timevar = "ind", direction = "wide")
            names(result) <- gsub(j, "exposure", names(result), fixed = TRUE)
            result$n <- nrow(na.omit(dat[, c(expo, i)]))

            # Add ANOVA p-values for main and interaction effects
            anova_ml <- anova(ml)
            result[, paste0("pval.aov.", rownames(anova_ml))] <- anova_ml$`Pr(>F)`

            return(result)
        })
        tmp <- do.call(rbind, tmp)
        return(tmp)
    })

    # Combine all results into one data frame
    res <- do.call(rbind, res)

    # Extract only essential information, including the sex:time interaction
    res <- res[, c("outcome", "n", grep("exposure|aov", names(res), value = TRUE), "pval.aov.sex:time")]

    return(res)
}

# Define parallelized group-specific function
mixed_anova_three_way_parallel_group <- function(dat, expo, outc, formel, group_var, covariates = NULL) {
    res <- foreach(i = seq_along(outc), .combine = rbind, .packages = c("lmerTest", "dplyr", "reshape2"),
                   .export = "mixed.anova.three_way_group") %dopar% {
                       mixed.anova.three_way_group(dat, expo, outc[i], formel, group_var, covariates)
                   }
    return(res)
}



# Example usage
# Define the data and parameters
# dat: your data frame
# expo: "t.factor" (time variable)
# outc: unique(prot.label$Assay) (list of proteins)
# formel: "+ (1|label)" (random effects and any additional terms)
# group_var: "group" (grouping variable)
# covariates: list of covariates, e.g., c("age", "sex", "total_fat_mass")




# Cluster temporal-profile plot for the Olink plasma proteome. Per-cluster
# protein trajectories keyed on OlinkID; the 5 lowest-q proteins are coloured
# (ggsci startrek), all others grey.
cluster_profile_olink <- function(res.exerome.linear, exerome.dat, prot.label, cl.label.exerome, cluster_number, ylim = c(-2.5, 6)) {
    mts <- res.exerome.linear %>%
        dplyr::filter(cluster == cluster_number) %>%
        dplyr::pull(OlinkID)
    filtered_data <- res.exerome.linear %>%
        dplyr::filter(cluster == cluster_number) %>%
        dplyr::select(OlinkID, Assay, fdr.aov)
    foo.m_list <- lapply(mts, function(m) {
        exerome.dat %>%
            dplyr::group_by(time_label) %>%
            dplyr::summarise(!!m := mean(!!sym(m)))
    })
    # Combine data into long format dataframe
    data <- bind_rows(lapply(seq_along(mts), function(i) {
        data.frame(
            OlinkID = rep(mts[i], nrow(foo.m_list[[i]])),
            Genes = filtered_data$Assay[filtered_data$OlinkID == mts[i]],
            time_label = foo.m_list[[i]]$time_label,
            zscore = foo.m_list[[i]][[2]],
            fdr.aov = filtered_data$fdr.aov[filtered_data$OlinkID == mts[i]], # Correcting the extraction of fdr.aov values
            cluster = cluster_number

        )
    }))
    # Arrange data based on fdr.aov in ascending order
    data <- data %>%
        arrange(fdr.aov)
    # Identify the lowest 5 unique entries
    lowest_5 <- data %>%
        distinct(Genes) %>%
        head(5)
    cl.vec <- c("Others" = "grey80", setNames(ggsci::pal_startrek()(6)[c(1:4, 6)][seq_len(nrow(lowest_5))], as.character(lowest_5$Genes)))
    # Create a new column to indicate whether each OlinkID is among the 5 lowest
    data <- data %>%
        group_by(Genes) %>%
        mutate(color_group = ifelse(Genes %in% lowest_5$Genes, as.character(Genes), "Others")) %>%
        ungroup()
    # Ensure that "Others" is placed at the bottom of the legend
    data$color_group <- factor(data$color_group, levels = c("Others", as.character(lowest_5$Genes)))
    # Create the plot using theming and aesthetics from generate_plot2
    p <- ggplot(data, aes(x = time_label, y = zscore, color = color_group, group = Genes)) +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = -Inf, ymax = Inf), fill = "grey95", color = NA, alpha = 0.7) +
        geom_line(data = data[data$color_group == "Others",], size = 0.3) + # Use size from generate_plot2
        geom_line(data = data[data$color_group != "Others",], size = 0.5) + # Use size from generate_plot2
        scale_color_manual(values = cl.vec) +
        geom_hline(yintercept = 0, lwd = 0.4, color = "black") +
        scale_x_continuous(breaks = c(-1, 0, 1, 3, 24), labels = c("-1", "0",  "1", "3", "24")) + # Match x scale
        labs(x = "Time [h]", y = "Z-score", title = paste0(cl.label.exerome$label[which(cl.label.exerome$cluster == cluster_number)], " - ", cl.label.exerome$cluster.paper[which(cl.label.exerome$cluster == cluster_number)])) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            text = element_text(size = 6),
            legend.key.size = unit(2, "mm"),
            panel.grid.minor = element_blank(),
            legend.position = c(1, 1),
            legend.justification = c("right", "top"),
            legend.background = element_blank(),
            plot.title = element_text(size = 6),
            axis.title = element_text(size = 6)
        ) +
        guides(color = guide_legend(title = NULL))  # Remove legend title

    if (!is.null(ylim)) p <- p + ggplot2::coord_cartesian(ylim = ylim)
    return(p)
}












# Cluster temporal-profile plot for the LC-MS/MS proteomes (saliva / urine /
# plasma LC-MS). Per-cluster protein trajectories keyed on Assay; the 5 lowest-q
# proteins are coloured (ggsci startrek), all others grey. Optional
# `save_data_path` writes the plotted long-format table to an .rds.
cluster_profile_lcms <- function(res.exerome.linear, exerome.dat, prot.label, cl.label.exerome, cluster_number, save_data_path = NULL, ylim = NULL) {
    mts <- res.exerome.linear %>%
        dplyr::filter(cluster == cluster_number) %>%
        dplyr::pull(Assay)
    filtered_data <- res.exerome.linear %>%
        dplyr::filter(cluster == cluster_number) %>%
        dplyr::select(Assay, fdr.aov)
    foo.m_list <- lapply(mts, function(m) {
        exerome.dat %>%
            dplyr::group_by(time_label) %>%
            dplyr::summarise(!!m := mean(!!sym(m), na.rm = TRUE))
    })
    data <- bind_rows(lapply(seq_along(mts), function(i) {
        data.frame(
            Assay = rep(mts[i], nrow(foo.m_list[[i]])),
            Genes = filtered_data$Assay[filtered_data$Assay == mts[i]],
            time_label = as.numeric(foo.m_list[[i]]$time_label),
            zscore = foo.m_list[[i]][[2]],
            fdr.aov = filtered_data$fdr.aov[filtered_data$Assay == mts[i]],
            cluster = cluster_number
        )
    })) %>%
        na.omit()

    if (!is.null(save_data_path)) {
        saveRDS(data, file = save_data_path)
        message("Dataframe saved to: ", save_data_path)
    }

    data <- data %>% arrange(fdr.aov)
    lowest_5 <- data %>% distinct(Genes) %>% head(5)
    cl.vec <- c("Others" = "grey80", setNames(ggsci::pal_startrek()(6)[c(1:4, 6)][seq_len(nrow(lowest_5))], as.character(lowest_5$Genes)))

    data <- data %>%
        group_by(Genes) %>%
        mutate(color_group = ifelse(Genes %in% lowest_5$Genes, as.character(Genes), "Others")) %>%
        ungroup()

    data$color_group <- factor(data$color_group, levels = c("Others", as.character(lowest_5$Genes)))

    p <- ggplot(data, aes(x = time_label, y = zscore, color = color_group, group = Genes)) +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = -Inf, ymax = Inf), fill = "grey95", color = NA, alpha = 0.7) +
        geom_line(data = data[data$color_group == "Others",], linewidth = 0.3) +
        geom_line(data = data[data$color_group != "Others",], linewidth = 0.5) +
        scale_color_manual(values = cl.vec) +
        geom_hline(yintercept = 0, lwd = 0.4, color = "black") +
        scale_x_continuous(breaks = c(-1, 0, 1, 3, 24), labels = c("-1", "0",  "1", "3", "24")) +
        labs(x = "Time [h]", y = "Z-score", title = paste0(cl.label.exerome$label[which(cl.label.exerome$cluster == cluster_number)], " - ", cl.label.exerome$cluster.paper[which(cl.label.exerome$cluster == cluster_number)])) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            text = element_text(size = 6),
            legend.key.size = unit(2, "mm"),
            panel.grid.minor = element_blank(),
            legend.position = c(1, 1),
            legend.justification = c("right", "top"),
            legend.background = element_blank(),
            plot.title = element_text(size = 6),
            axis.title = element_text(size = 6)
        ) +
        guides(color = guide_legend(title = NULL))

    if (!is.null(ylim)) p <- p + ggplot2::coord_cartesian(ylim = ylim)
    return(p)
}


# exerome generate plot

metabolic_markers <- function(res.exerome.linear, exerome.dat, prot.label, cl.label.exerome, cluster_number, ylim = c(-2.5, 6)) {
    mts <- res.exerome.linear %>%
        dplyr::filter(cluster == cluster_number) %>%
        dplyr::pull(Assay)

    filtered_data <- res.exerome.linear %>%
        dplyr::filter(cluster == cluster_number) %>%
        dplyr::select(Assay, pval.aov.t.factor)

    foo.m_list <- lapply(mts, function(m) {
        exerome.dat %>%
            dplyr::group_by(time_label) %>%
            dplyr::summarise(!!m := mean(!!sym(m), na.rm = TRUE))  # Handle NA values
    })

    # Combine data into long format dataframe
    data <- bind_rows(lapply(seq_along(mts), function(i) {
        data.frame(
            Assay = rep(mts[i], nrow(foo.m_list[[i]])),
            Genes = filtered_data$Assay[filtered_data$Assay == mts[i]],
            time_label = as.numeric(foo.m_list[[i]]$time_label),  # Ensure time_label is numeric
            zscore = foo.m_list[[i]][[2]],
            pval.aov.t.factor = filtered_data$pval.aov.t.factor[filtered_data$Assay == mts[i]],  # Correcting the extraction of pval.aov.t.factor values
            cluster = cluster_number
        )
    })) %>%
        na.omit()  # Remove rows with NA values

    # Create a new column for color values and use original names for assays
    data <- data %>%
        mutate(
            color_group = ifelse(pval.aov.t.factor < 0.05, Genes, "non-significant")
        )

    # Generate a range of grey colors for non-significant assays
    non_sig_genes <- unique(data$Genes[data$color_group == "non-significant"])
    non_sig_colors <- setNames(colorRampPalette(c("grey90", "grey80"))(length(non_sig_genes)), non_sig_genes)

    # Set color values for the plot
    sig_genes <- unique(data$color_group[data$color_group != "non-significant"]); sig_colors <- setNames(grDevices::colorRampPalette(ggsci::pal_startrek()(6))(length(sig_genes)), sig_genes)
    cl.vec <- c(sig_colors, non_sig_colors)

    # Sort data by pval.aov.t.factor to determine the order of the legend
    data <- data %>%
        arrange(pval.aov.t.factor)

    # Set the order of color_group levels to ensure proper legend ordering
    data <- data %>%
        mutate(color_group = factor(Genes, levels = unique(Genes)))

    # Create the plot using theming and aesthetics from generate_plot2
    p <- ggplot(data, aes(x = time_label, y = zscore, color = color_group, group = Genes)) +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = -Inf, ymax = Inf), fill = "grey95", color = NA, alpha = 0.7) +
        geom_line(data = data, linewidth = 0.5) +  # Use linewidth instead of size
        scale_color_manual(values = cl.vec) +
        geom_hline(yintercept = 0, lwd = 0.4, color = "black") +
        scale_x_continuous(breaks = c(-1, 0, 1, 3, 24), labels = c("-1", "0",  "1", "3", "24")) +  # Match x scale
        labs(x = "Time [h]", y = "Z-score", title = cl.label.exerome$label[which(cl.label.exerome$cluster == cluster_number)]) +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            text = element_text(size = 6),
            legend.key.size = unit(2, "mm"),
            panel.grid.minor = element_blank(),
            legend.position = c(1, 1),
            legend.justification = c("right", "top"),
            legend.background = element_blank(),
            plot.title = element_text(size = 6),
            axis.title = element_text(size = 6)
        ) +
        guides(color = guide_legend(title = NULL))  # Remove legend title

    if (!is.null(ylim)) p <- p + ggplot2::coord_cartesian(ylim = ylim)
    return(p)
}

# this funciton is working for clusters with proteins with difficult NA values where lines arnt showing in plots
