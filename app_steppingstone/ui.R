
shinyUI(pageWithSidebar(
  headerPanel("Stepping stone process to model animal movement"),
  sidebarPanel(
    conditionalPanel(condition="input.conditionedPanels==1",
                     numericInput("nc", "Number of Cells", 200, 50, 500),
                     h4("Attraction to home-range center"),
                     helpText("Squared distance to the home-range center."),
                     numericInput("cent_x", "Home-range center (x)", 100, 50, 500),
                     numericInput("cent_y", "Home-range center (y)", 100, 50, 500),
                     h4("Habitat"),
                     helpText("Paramaeters for Gauss Random Fields"),
                     sliderInput("sig", "Partial sill (sigma^2)", 0, 100, 1),
                     sliderInput("phi", "Range (phi)", 0, 100, 10)
    ),
    conditionalPanel(condition="input.conditionedPanels==2",
                     sliderInput("n", "No steps", 1e1, 1e6, 1e5),
                     sliderInput("rari", "Rarify by", 1, 1e3, 10),
                     sliderInput("burnin", "Burnin", 0, 1e4, 0),
                     sliderInput("pm", "Probability of moving", 0.0001, 0.9999, 0.5),
                     sliderInput("omg_dist", "Omega (dist home-range center)", -0.5, 0, -0.05, step = 0.0001),
                     sliderInput("omg_hab", "Omega (habitat)", -10, 10, 1, step = 0.01),
                     sliderInput("kappa", "Kappa (VonMises)", 0, 20, 1, step = 0.1),
                     selectInput("boundary", "Boundary condition", choices = c("Wrapped" = "wrap",
                                                                               "Stop at boundary" = "stop",
                                                                               "Reflective" = "reflective")),
                     numericInput("st_x", "Start (x)", 50),
                     numericInput("st_y", "Start (y)", 50)
    )
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Covariates", value=1,
               plotOutput("dist_hrc"),
               plotOutput("hab")
      ),
      tabPanel("Movement", value=2,
               h2("Adjust Plot"),
               selectInput("background", "Background", c("none", "dist to home-range center" = "dist_hrc", "habitat" = "hab")),
               sliderInput("trans", "Transparency", 0, 1, 0.5, 0.01),
               checkboxInput("start", "Show start"),
               plotOutput("w")
      ),
       id = "conditionedPanels"
    )
  )
))
