
# load rds table final_dt

final_dt[, class_plot := fifelse(
    !is.na(classification),
    classification,
    fifelse(
        !is.na(variant.id_mel) & is.na(variant.id_sim),
        "mel_only",
        fifelse(
            is.na(variant.id_mel) & !is.na(variant.id_sim),
            "sim_only",
            NA_character_
        )
    )
)]

table(final_dt$class_plot, useNA = "ifany")

library(ggplot2)

class_counts <- final_dt[
    !is.na(class_plot),
    .N,
    by = class_plot
]

class_counts[, class_plot := factor(
    class_plot,
    levels = c(LETTERS[1:25], "mel_only", "sim_only")
)]

ggplot(class_counts, aes(x = class_plot, y = N)) +
    geom_col() +
    geom_text(
        aes(
            label = N,
            y = pmin(N, 70000)
        ),
        angle = 45,
        hjust = 0,
        vjust = -1, 
        size = 2
    ) +
    coord_cartesian(ylim = c(0, 80000)) +
    labs(
        x = "Classification",
        y = "Number of SNPs",
        title = "Distribution of SNPs by classification"
    ) +
    theme_classic()