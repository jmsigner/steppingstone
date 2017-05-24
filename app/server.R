library(shiny)
library(steppingstone)
library(geoR)
library(raster)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {

  r <- reactive({
    r <- raster(xmn = 0, xmx = input$nc, ymn = 0, ymx = input$nc, res = 1)
    # dist hr-center
    r <- distanceFromPoints(r, cbind(input$cent_x, input$cent_y))
    (r^2)/10
  })

  h <- reactive({
    raster(geoR::grf(ncell(r()), grid="reg", nx=input$nc, ny=input$nc, xlims=c(0, input$nc), ylims=c(0, input$nc),
                    cov.pars=c(input$sig, input$phi)))

  })

  # simulate movement
  w <- reactive({
    rsc <- stack(r(), h())
    w <- steppingstone(alpha = pm2alpha(input$pm),
                       omegas = c(input$omg_dist, input$omg_hab), rsc,
                       kappa = input$kappa,
                       xy0 = c(input$st_x, input$st_y),
                       n = input$n,
                       burnin = input$burnin,
                       boundary = input$boundary,
                       rarify_by = input$rari)
  })

  output$dist_hrc <- renderPlot(
    plot(r())
  )

  output$hab <- renderPlot(
    plot(h())
  )

  output$w <- renderPlot({
    if (input$background == "none") {
      plot(-10, -10, xlim = c(0, input$nc), ylim = c(0, input$nc), xlab = "", ylab = "", asp = 1)
    } else if (input$background == "hab") {
      plot(h())
    } else {
      plot(r())
    }
    points(w()$xy, pch = 20, col = adjustcolor("black", alpha = input$trans))

    if (input$start) {
      points(w()$xy[1, 1], w()$xy[1, 2], pch = 20, cex = 3, col = "red")
    }
  })

})
