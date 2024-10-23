# How many people was infected ----
# Histogram of how many people infected do we have in each wave
t <- waves_data |>
  select("eid", ends_with("_antibody_test_result")) |>
  filter(across(-"eid", ~!is.na(.))) 

# People infected at each wave:
library(tidyr)
people_infected_in_each_wave <- t |> 
  select(-c("eid")) |>
  summarise(across(everything(), ~sum(.))) |>
  pivot_longer(cols = everything()) |>
  mutate(name = gsub("_antibody_test_result","",name))

# How many people infected at least one time
w0 <- t |>
  select("eid", "w0_antibody_test_result") |>
  filter(w0_antibody_test_result == 1)
w1 <- t |>
  select("eid", "w1_antibody_test_result") |>
  filter(w1_antibody_test_result == 1) |>
  anti_join(w0 |> select("eid"))
w2 <- t |>
  select("eid", "w2_antibody_test_result") |>
  filter(w2_antibody_test_result == 1) |>
  anti_join(w0 |> select("eid")) |>
  anti_join(w1 |> select("eid"))
w3 <- t |>
  select("eid", "w3_antibody_test_result") |>
  filter(w3_antibody_test_result == 1) |>
  anti_join(w0 |> select("eid")) |>
  anti_join(w1 |> select("eid")) |>
  anti_join(w2 |> select("eid"))
w4 <- t |>
  select("eid", "w4_antibody_test_result") |>
  filter(w4_antibody_test_result == 1) |>
  anti_join(w0 |> select("eid"))|>
  anti_join(w1 |> select("eid"))|>
  anti_join(w2 |> select("eid"))|>
  anti_join(w3 |> select("eid"))
w5 <- t |>
  select("eid", "w5_antibody_test_result") |>
  filter(w5_antibody_test_result == 1)|>
  anti_join(w0 |> select("eid"))|>
  anti_join(w1 |> select("eid"))|>
  anti_join(w2 |> select("eid"))|>
  anti_join(w3 |> select("eid"))|>
  anti_join(w4 |> select("eid"))

new_people_infected_in_each_wave <- tibble(
  w0 = w0 |> nrow(),
  w1 = w1 |> nrow(), 
  w2 = w2 |> nrow(),
  w3 = w3 |> nrow(),
  w4 = w4 |> nrow(),
  w5 = w5 |> nrow(),
)




