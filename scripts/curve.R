library(lmtp)
library(ggplot2)
library(purrr)
library(dplyr)

odtr <- readRDS("data/drv/survival-combined-odtr.rds")
simple <- readRDS("data/drv/survival-combined-dynamic.rds")

combined <- readRDS("data/src/clean•weeks•with•relapse•wide•010422.rds")

bup <- filter(combined, medicine == "bup")
Y <- glue::glue("wk{4:12}.relapse_this_week")

obs <- tibble::tibble(
    time = 4:12, 
    strategy = "Observed", 
    estimate = map_dbl(Y, \(x) mean(bup[[x]]))
)

tidy(odtr) |> 
    mutate(estimate = 1 - estimate, 
           conf.low = estimate + qnorm(0.025)*std.error, 
           conf.high = estimate - qnorm(0.025)*std.error)

tidy(simple) |> 
    mutate(estimate = 1 - estimate, 
           conf.low = estimate + qnorm(0.025)*std.error, 
           conf.high = estimate - qnorm(0.025)*std.error)

# odtr <- 
#     map_dfr(odtr[4:12], tidy) |> 
#     mutate(estimate = if_else(row_number() != 1, 1 - estimate, estimate), 
#            conf.low = estimate + qnorm(0.025)*std.error, 
#            conf.high = estimate - qnorm(0.025)*std.error) |> 
#     mutate(time = 4:12, .before = 2) |> 
#     mutate(strategy = "ODTR", .before = 3)
# 
# simple <- 
#     map_dfr(simple[3:11], tidy) |> 
#     mutate(estimate = if_else(row_number() != 1, 1 - estimate, estimate), 
#            conf.low = estimate + qnorm(0.025)*std.error, 
#            conf.high = estimate - qnorm(0.025)*std.error) |> 
#     mutate(time = 4:12, .before = 2) |> 
#     mutate(strategy = "Simple", .before = 3)

ragg::agg_png("figures/combined-survival-curve.png", width = 8, height = 4.5, units = "cm", res = 400)
bind_rows(obs, simple, odtr) |> 
ggplot(aes(x = time, y = estimate, color = strategy)) + 
    geom_step(size = 0.2) + 
    geom_point(size = 0.2, aes(shape = strategy, color = strategy)) + 
    scale_x_continuous(
        breaks = 4:12, 
        labels = c("Wk. 4", 5:12), 
        limits = c(3.75, 12.25), 
        expand = c(0, .2)
    ) + 
    scale_color_manual(values = c('#004488', '#DDAA33', '#BB5566')) + #, '#000000')) 
    labs(
        x = "", 
        y = "Relapse risk", 
        shape = "", 
        color = ""
    ) + 
    theme_classic(
        base_size = 4, 
        base_line_size = 0.2, 
        base_rect_size = 0.2
    ) + 
    theme(
        legend.text = element_text(size = 3), 
        legend.key.size = unit(0.25, "cm")
    )  + 
    guides(guide_legend(override.aes = list(size = 0.5)))
dev.off()
