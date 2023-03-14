
#' Title
#'
#' @return
#' @export
#'
#' @examples
#' @import caribouMetrics shiny dplyr shinydashboard
#'
demographicProjectionApp <- function() {
  rjags::load.module("glm")

  # TO DO: make recruitment adjustment (equation 2 of Eacker et al) an option. Note this requires adjusting all calls to popGrowthJohnson, as well as jags code.

  ##########
  # Get full set of sims for comparison
  # TO DO: need an option for force update of this cache - needed when something changes in caribouMetrics, which presumably won't happen often, but may happen sometimes.
  # Include this option with a note that it is something to try if one is getting a mismatch between local and national results given minimal monitoring data.
  simBigAdjust <- getSimsNational(adjustR = T) # If called with default parameters, use saved object to speed things up.
  simBigNoAdjust <- getSimsNational(adjustR = F) # If called with default parameters, use saved object to speed things up.

  # get default values to use in initializing parameters in ui
  # Disturbance scenario parameters
  # P: num obs years
  # J: num projection years
  # aS: % increase in anthropogenic disturbance per year in observation period
  # aSf: % increase in anthropogenic disturbanc per year in projection period
  # iA: initial anthro disturbance
  # iF: initial fire - for now, no change in fire over time.
  # Or, if disturbance csv scenario provided, ignore all these setting and get values from file.

  # True population parameters
  # N0: intial population size
  # rQ: recruitment quantile
  # sQ: survival quantile
  # rS: disturbance-recruitment slope multiplier
  # sS: disturbance-survival slope multiplier
  # adjustR: Adjust R to account for delayed age at first reproduction (DeCesare et al. 2012; Eacker et al. 2019) or not.
  # In theory, could expose any parameter from caribouMetrics projection tool here. Leave simple for now.

  # Obs model parameters
  # st: target number of collars per year
  # ri: renewal interval

  scn_defaults <- c(eval(formals(fillDefaults)$defList),
    curYear = eval(formals(fillDefaults)$curYear)
  )

  # Observation model parameters
  # cowsPerYear: number of cows in aerial surval for calf:cow ratio each year.
  # startsPerYear: # collars deployed each year
  # collarNumYears: years before each collar falls off
  # collarOnTime: month that collars are deployed
  # collarOffTime: month that collars fall off

  # defaults set to be uninformative
  obs_defaults <- list(
    cmult = 0, startsPerYear = 1, renewalInterval = 1, collarNumYears = 6,
    collarOnTime = 1, collarOffTime = 12
  )

  # Priors - see getPriors() for details & defaults
  # Allow possibility of a different parameter table instead of caribouMetrics::popGrowthTableJohnsonECCC. But don't allow different covariates. For now, recruitment model must have anthro and fire_excl_anthro covariates, and survival model must have anthro covariate. Intercept and slope parameters must be specified for each.
  # bse: anthro slope survival uncertainty multiplier. 1 - 10
  # bre: anthro slope recruitment uncertainty multiplier. 1 - 10
  # lse: survival intercept uncertainty multiplier. 1 - 10
  # lre: recruitment intercept uncertainty multiplier. 1 - 10
  # sse: interannual coefficient of variation for survival. 0-1. See popGrowthJohnson() and functions therein for details
  # ssv: uncertainty about interannual variation in survival. 0-1
  # sre: interannual coefficient of variation for recruitment. 0-1
  # srv: uncertainty about interannual variation in recruitment. 0-1.
  prior_defaults <- eval(formals(getPriors)$expectMods)

  # JAGS params
  jags_defaults <- eval(formals(runRMModel)[c("Niter", "Nchains", "Nthin", "Nburn")])

  # add JavaScript to add an id to the <section> tag
  # so we can overlay waiter on top of it
  add_id_to_section <- "
$( document ).ready(function() {
  var section = document.getElementsByClassName('content');
  section[0].setAttribute('id', 'waiter-content');
});"

  # UI  -------------------------------------
  ui <- dashboardPage(
    # Application title
    dashboardHeader(title = "Woodland Caribou Demographic Estimates"),

    # Sidebar --------------------------------------------------------------------
    dashboardSidebar(
      width = 400,
      h2("Define simulation inputs"),
      shinyjs::useShinyjs(),
      sidebarMenu(
        id = "sideBar",
        # disturbance --------------------------------
        menuItem(
          text = "Scenario",
          selectInput("scn_source",
            label = "Create disturbance scenario from...",
            choices = c("a csv file" = "file", "a simulation" = "sim"),
            selected = "sim"
          ),
          numericInput("curYear",
            label = "Last year of observations",
            value = scn_defaults$curYear
          ),
          shinyFiles::shinyFilesButton("scn_file", "Select File", "Please select a file",
            multiple = FALSE
          ),
          div(
            id = "nYears",
            numericInput("nYearObs",
              label = "Number of years of observations",
              value = scn_defaults$P, min = 1
            ),
            numericInput("nYearProj",
              label = "Number of years of projections",
              value = scn_defaults$J, min = 0
            )
          ),
          div(
            id = "dist_sim",
            numericInput("iA",
              label = "Initial % anthropogenic disturbance",
              value = scn_defaults$iA, min = 0
            ),
            numericInput("aS",
              label = "% increase in anthropogenic disturbance per year in observation period",
              value = scn_defaults$aS, min = 0
            ),
            numericInput("aSf",
              label = "% increase in anthropogenic disturbance per year in projection period",
              value = scn_defaults$aSf, min = 0
            ),
            numericInput("iF",
              label = "Initial % natural disturbance",
              value = scn_defaults$iF, min = 0
            )
          )
        ),
        # True pop ---------------------------
        menuItem(
          "True population parameters",
          checkboxInput("adjustR", "Adjust R to account for delayed age at first reproduction", value = scn_defaults$adjustR),
          numericInput("N0",
            label = "Initial population size",
            value = scn_defaults$N0, min = 0
          ),
          sliderInput(
            inputId = "rQ", label = "Recruitment quantile",
            value = scn_defaults$rQ, min = 0.025, max = 0.975
          ),
          sliderInput(
            inputId = "sQ", label = "Survival quantile",
            value = scn_defaults$sQ, min = 0.025, max = 0.975
          ),
          sliderInput(
            inputId = "rS", label = "Multiplier for effect of disturbance on recruitment",
            value = scn_defaults$rS, min = 0, max = 5
          ),
          sliderInput(
            inputId = "sS", label = "Multiplier for effect of disturbance on survival",
            value = scn_defaults$sS, min = 0, max = 5
          )
        ),
        # Obs data ---------------------------------
        menuItem(
          "Observation model parameters",
          numericInput("startsPerYear",
            label = "Target number of collars",
            value = obs_defaults$startsPerYear, min = 1
          ),
          numericInput("renewalInterval",
            label = "Number of years between collar deployments",
            value = obs_defaults$renewalInterval, min = 1
          ),
          numericInput("cmult",
            label = "Number of cows per collared cow in aerial surveys for calf:cow ratio each year",
            value = obs_defaults$cmult, min = 0
          ),
          numericInput("collarNumYears",
            label = "Number of years until collar falls off",
            value = obs_defaults$collarNumYears, min = 1
          ),
          numericInput("collarOnTime",
            label = "Month that collars are deployed",
            value = obs_defaults$collarOnTime, min = 1, max = 12
          ),
          numericInput("collarOffTime",
            label = "Month that collars fall off",
            value = obs_defaults$collarOffTime, min = 1, max = 12
          )
        ),
        # Priors ------------------------------------------------
        menuItem(
          "Model priors",
          checkboxInput("redoSimsNational", "Update cached national simulations.", value = 0),
          selectInput("nat_model",
            label = "National version model to use",
            choices = c("default", "custom"),
            selected = "default"
          ),
          div(
            id = "custModel",
            h4(HTML("To use a different version of the national<br/>model select the model numbers below")),
            selectInput("recruitmentModelNumber",
              label = "Recruitment model number",
              choices = caribouMetrics::popGrowthTableJohnsonECCC %>%
                filter(responseVariable == "recruitment") %>%
                group_by(ModelNumber) %>%
                filter(all(Coefficient %in% c(
                  "Intercept", "Precision",
                  "Anthro", "fire_excl_anthro"
                ))) %>%
                pull(ModelNumber) %>% unique(),
              selected = "M4"
            ),
            selectInput("survivalModelNumber",
              label = "Female survival model number",
              choices = caribouMetrics::popGrowthTableJohnsonECCC %>%
                filter(responseVariable == "femaleSurvival") %>%
                group_by(ModelNumber) %>%
                filter(all(Coefficient %in% c(
                  "Intercept", "Precision",
                  "Anthro", "fire_excl_anthro"
                ))) %>%
                pull(ModelNumber) %>% unique(),
              selected = "M1"
            ),
            h4(HTML("To use a custom table of coefficients<br/>select a csv file")),
            shinyFiles::shinyFilesButton("popGrowth_file", "Select File", "Please select a file",
              multiple = FALSE
            )
          ),
          sliderInput(
            inputId = "bse",
            label = "Multiplier for uncertainty of effect of disturbance on survival",
            value = prior_defaults$bse, min = 1, max = 10
          ),
          sliderInput(
            inputId = "bre",
            label = "Multiplier for uncertainty of effect of disturbance on recruitment",
            value = prior_defaults$bre, min = 1, max = 10
          ),
          sliderInput(
            inputId = "lse",
            label = "Multiplier for uncertainty of survival intercept",
            value = prior_defaults$lse, min = 1, max = 10
          ),
          sliderInput(
            inputId = "lre",
            label = "Multiplier for uncertainty of recruitment intercept",
            value = prior_defaults$lre, min = 1, max = 10
          ),
          sliderInput(
            inputId = "sse",
            label = "Interannual coefficient of variation for survival",
            value = prior_defaults$sse, min = 0, max = 1
          ),
          sliderInput(
            inputId = "ssv",
            label = "Uncertainty around interannual coefficient of variation for survival",
            value = prior_defaults$ssv, min = 0, max = 1
          ),
          sliderInput(
            inputId = "sre",
            label = "Interannual coefficient of variation for recruitment",
            value = prior_defaults$sre, min = 0, max = 1
          ),
          sliderInput(
            inputId = "srv",
            label = "Uncertainty around interannual coefficient of variation for recruitment",
            value = prior_defaults$srv, min = 0, max = 1
          )
        ),

        # JAGS params ---------------------------
        menuItem(
          "Baysian model parameters",
          checkboxInput("getKSDists", "Calculate Kolmogorov-Smirnov Distances", value = 0),
          selectInput("survAnalysisMethod",
            label = "Survival analysis method to use",
            choices = c("KaplanMeier", "Exponential"),
            selected = "KaplanMeier"
          ),
          sliderInput(
            inputId = "Nchains", label = "Number of chains",
            value = jags_defaults$Nchains, min = 1, max = 5
          ),
          sliderInput(
            inputId = "Niter", label = "Number of iterations",
            value = jags_defaults$Niter, min = 1, max = 50000, step = 1000
          ),
          sliderInput(
            inputId = "Nburn", label = "Length of burn-in",
            value = jags_defaults$Nburn, min = 1, max = 20000, step = 1000
          ),
          sliderInput(
            inputId = "Nthin", label = "Thinning rate",
            value = jags_defaults$Nthin, min = 1, max = 10
          )
        )
        # ------------------------
      ),
      h3("Run model/Update data"),
      actionButton("Run.model", "Run model", icon("paper-plane"),
        style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
      ),
      actionButton("mcmc.output", "Open MCMC diagnostic plots"),
      shinyFiles::shinySaveButton("saveOutput",
        title = "Save results to csv or rds",
        label = "Save results to csv or rds",
        filetype = c("csv", "rds")
      )
      # This .rmd file is missing
      # radioButtons('format', 'Document format', c('HTML', 'PDF', 'Word'),
      #              selected="HTML",  inline = TRUE),
      # downloadButton("downloadOutput", label = "Download output"),
    ),
    # Body #---------------------------------------------------------------------
    dashboardBody(
      # import our custom JavaScript
      tags$head(
        tags$script(add_id_to_section)
      ),
      waiter::useWaiter(),
      navbarPage(
        id = "bodyTabs",
        title = "",
        tabPanel(
          "Instructions",
          includeMarkdown(system.file("app/Instructions.md",
            package = "BayesianCaribouDemographicProjection"
          ))
        ),
        tabPanel(
          "Graphical summary",
          tabsetPanel(
            id = "graphPanel",
            tabPanel("Disturbance", plotOutput("plot6")),
            tabPanel("Recruitment", plotOutput("plot2")),
            tabPanel("Adult female survival", plotOutput("plot1")),
            tabPanel("Population growth rate", plotOutput("plot4")),
            tabPanel("Female population size", plotOutput("plot5")),
            tabPanel("Female recruitment", plotOutput("plot3")),
            tabPanel("Recruitment Kolmogorov-Smirnov Distance", plotOutput("plot8")),
            tabPanel("Adult female survival Kolmogorov-Smirnov Distance", plotOutput("plot7")),
            tabPanel("Population growth rate Kolmogorov-Smirnov Distance", plotOutput("plot9"))
          ),
        ),
        tabPanel(
          "Tabular summary",
          tabsetPanel(
            id = "tablePanel",
            tabPanel("Summary", dataTableOutput("table")),
            tabPanel("Adult female survival", tableOutput("table2")),
            tabPanel("Recruitment", tableOutput("table3")),
            tabPanel("Female recruitment", tableOutput("table4")),
            tabPanel("Population growth rate", tableOutput("table5")),
            tabPanel("Female population size", tableOutput("table7")),
            tabPanel("JAGS output", dataTableOutput("table6"))
          ),
        )
      )
    )
  )



  # Define server logic
  server <- function(input, output, session) {

    # show appropriate inputs based on file or sim scenario source
    observe({
      switch(input$scn_source,
        sim = {
          shinyjs::hide("scn_file")
          shinyjs::show("nYears")
          shinyjs::show("dist_sim")
        },
        file = {
          shinyjs::show("scn_file")
          shinyjs::hide("dist_sim")
          shinyjs::hide("nYears")
        }
      )

      switch(input$nat_model,
        default = {
          shinyjs::hide("custModel")
        },
        custom = {
          shinyjs::show("custModel")
        }
      )
    })

    volumes <- eventReactive(c(input$scn_file, input$saveOutput, input$popGrowth_file),
      ignoreInit = TRUE,
      {
        message(
          "Accessing local file structure.",
          "\nIf this takes a long time it is because you connected and then ",
          "disconnected from VPN. Either reconnect or restart your computer to ",
          "speed up"
        )

        # Set up file selection
        # Note this time out is because when I disconnected from VPN it
        # made the getVolumes function hang forever because it was looking for
        # drives that were no longer connected. Now it will give an error
        timeout <- R.utils::withTimeout(
          {
            volumes <- c(
              wd = getwd(),
              getVolumes()()
            )
          },
          timeout = 200,
          onTimeout = "silent"
        )

        if (is.null(timeout)) {
          stop("The app is unable to access your files because you were connected",
            " to the VPN and then disconnected. To fix this either reconnect to",
            " the VPN or restart your computer and use the app with out connecting",
            " to VPN. See issue https://github.com/see24/ccviR/issues/36 for more ",
            "information",
            call. = FALSE
          )
        }

        return(volumes)
      }
    )

    observeEvent(c(input$scn_file, input$saveOutput, input$popGrowth_file), {
      shinyFiles::shinyFileChoose(input, "scn_file",
        roots = volumes(), session = session,
        filetypes = "csv"
      )

      shinyFiles::shinyFileChoose(input, "popGrowth_file",
        roots = volumes(), session = session,
        filetypes = "csv"
      )

      shinyFiles::shinyFileSave(input, "saveOutput",
        roots = volumes(), session = session,
        filetypes = c("csv", "rds")
      )
    })

    scn_df <- eventReactive(input$scn_file, {
      # if(!is.integer(input$scn_file)){
      df <- read.csv(parseFilePaths(volumes, input$scn_file)$datapath)
      # check colnames match expected for caribouMetrics
      missed_nms <- setdiff(c("Anthro", "fire_excl_anthro", "Year"), colnames(df))
      if (length(missed_nms) > 0) {
        stop("The scenario file loaded is missing the columns ",
          paste(missed_nms, collapse = ", "),
          call. = FALSE
        )
      }

      return(df)
    })

    popGrow_df <- eventReactive(input$popGrowth_file, {
      if (!is.integer(input$popGrowth_file)) {
        df <- read.csv(parseFilePaths(volumes, input$popGrowth_file)$datapath,
          blank.lines.skip = TRUE
        )
        # check table matches expected for caribouMetrics
        # TODO: confirm that a table that passes this works
        df <- testPopGrowthTable(df)

        return(df)
      }
    })
    observe(print(popGrow_df()))

    observeEvent(input$scn_file, {
      if (!is.integer(input$scn_file)) {
        show("nYears")
        # set the inputs to these values based on the input data and enforce a max
        updateNumericInput(
          inputId = "nYearObs",
          value = input$curYear - min(scn_df()$Year) + 1,
          max = input$curYear - min(scn_df()$Year) + 1
        )

        updateNumericInput(
          inputId = "nYearProj",
          value = max(scn_df()$Year) - input$curYear,
          max = max(scn_df()$Year) - input$curYear
        )
      }
    })

    waiter <- waiter::Waiter$new(id = c("waiter-content"))

    observeEvent(input$Run.model, updateNavbarPage(
      inputId = "bodyTabs",
      selected = "Graphical summary"
    ))


    dataInput <- eventReactive(input$Run.model, {

      # TO DO: consider making this happen when update box is checked?
      if (input$redoSimsNational) {
        simBigAdjust <- getSimsNational(adjustR = T, forceUpdate = T, fire_excl_anthro = input$iF)
        simBigNoAdjust <- getSimsNational(adjustR = F, forceUpdate = T, fire_excl_anthro = input$iF)
      }

      waiter$show()
      waiter$update(html = tagList(
        waiter::spin_1(),
        h4("Running model...")
      ))
      on.exit(waiter$hide())

      startYear <- input$curYear - input$nYearObs + 1

      endYear <- input$curYear + input$nYearProj

      scns <- expand.grid(
        P = input$nYearObs, J = input$nYearProj,
        aS = input$aS, aSf = input$aSf, rS = input$rS,
        sS = input$sS, iA = input$iA, iF = input$iF,
        rQ = input$rQ, sQ = input$sQ, N0 = input$N0,
        iYr = startYear, adjustR = input$adjustR, cmult = input$cmult, st = input$startsPerYear, ri = input$renewalInterval
      )


      scns <- fillDefaults(scns)

      # create cowCounts and freqStartsByYear from inputs
      cowCounts <- data.frame(
        Year = startYear:input$curYear,
        Count = input$startsPerYear * input$cmult,
        Class = "cow"
      )
      freqStartsByYear <- data.frame(
        Year = startYear:input$curYear,
        numStarts = input$startsPerYear
      )
      if (is.integer(input$scn_file)) {
        scn_df2 <- NULL
      } else {
        scn_df2 <- scn_df()
      }

      if (input$nat_model == "default") {
        popGrow_df2 <- caribouMetrics::popGrowthTableJohnsonECCC
      } else {
        if (is.integer(input$popGrowth_file)) {
          popGrow_df2 <- caribouMetrics::popGrowthTableJohnsonECCC
        } else {
          popGrow_df2 <- popGrow_df()
        }
      }

      oo <- simulateObservations(scns,
        distScen = scn_df2, cowCounts = cowCounts,
        freqStartsByYear = freqStartsByYear,
        collarNumYears = input$collarNumYears,
        collarOffTime = input$collarOffTime,
        collarOnTime = input$collarOnTime,
        populationGrowthTable = popGrow_df2,
        survivalModelNumber = input$survivalModelNumber,
        recruitmentModelNumber = input$recruitmentModelNumber
      )

      betaPriors <- getPriors(
        modifiers = isolate(reactiveValuesToList(input)),
        popGrowthTable = popGrow_df2
      )

      out <- runRMModel(
        survData = oo$simSurvObs, ageRatio.herd = oo$ageRatioOut,
        disturbance = oo$simDisturbance, betaPriors = betaPriors,
        startYear = startYear, endYear = endYear,
        inpFixed = input
      )
      out$oo <- oo
      out$startYear <- startYear
      out$endYear <- endYear
      out
    })


    # TABLES #######
    # TODO: there is a simpler way to do this reactivity with req()
    dataInput1 <- eventReactive(input$Run.model, {
      out <- dataInput()
      if (input$adjustR) {
        simBig <- simBigAdjust
      } else {
        simBig <- simBigNoAdjust
      }
      return(getOutputTables(
        result = out$result, startYear = out$startYear, endYear = out$endYear,
        survInput = out$survInput, oo = out$oo, simBig = simBig, getKSDists = input$getKSDists
      ))
    })

    output$table <- renderDataTable({
      dataInput1()$rr.summary.all %>%
        # remove params as they are included in label
        select(-c(P, J, aS, aSf, rS, sS, iA, iF, rQ, sQ, N0, iYr, st, ID)) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2)))
    })

    dataInput2 <- eventReactive(input$Run.model, {
      df <- getSumStats("S.annual.KM", dataInput()$result, dataInput()$startYear, dataInput()$endYear)

      names(df)[3] <- "Survival"

      df
    })

    output$table2 <- renderTable({
      dataInput2()
    })

    dataInput3 <- eventReactive(input$Run.model, {
      getSumStats("R", dataInput()$result, dataInput()$startYear, dataInput()$endYear)
    })

    output$table3 <- renderTable({
      dataInput3()
    })

    dataInput4 <- eventReactive(input$Run.model, {
      getSumStats("Rfemale", dataInput()$result, dataInput()$startYear, dataInput()$endYear)
    })

    output$table4 <- renderTable({
      dataInput4()
    })

    dataInput5 <- eventReactive(input$Run.model, {
      getSumStats("pop.growth", dataInput()$result, dataInput()$startYear, dataInput()$endYear)
    })

    output$table5 <- renderTable({
      dataInput5()
    })

    dataInput7 <- eventReactive(input$Run.model, {
      getSumStats("fpop.size", dataInput()$result, dataInput()$startYear, dataInput()$endYear)
    })

    output$table7 <- renderTable({
      dataInput7()
    })

    output$table6 <- renderDataTable({
      columnNAMES <- colnames(dataInput()$result$BUGSoutput$summary)

      x <- data.frame(
        row.names(dataInput()$result$BUGSoutput$summary),
        dataInput()$result$BUGSoutput$summary
      )

      names(x) <- c("Parameters", columnNAMES)
      x
    })

    observeEvent(input$mcmc.output, {
      req(dataInput())
      p <- mcmcplot(dataInput()$result)
      print(p)
    })


    # PLOTS #######

    # Adult female survival
    output$plot1 <- renderPlot({
      scResults <- dataInput1()
      plotRes(scResults$rr.summary.all, "Adult female survival",
        obs = scResults$obs.all,
        lowBound = 0.6, simRange = scResults$sim.all
      )
    })

    # Recruitment
    output$plot2 <- renderPlot({
      scResults <- dataInput1()
      plotRes(scResults$rr.summary.all, "Recruitment",
        obs = scResults$obs.all,
        lowBound = 0, simRange = scResults$sim.all
      )
    })

    # Female-only Recruitment
    output$plot3 <- renderPlot({
      plotRes(dataInput1()$rr.summary.all, "Female-only recruitment")
    })

    # lambda
    output$plot4 <- renderPlot({
      scResults <- dataInput1()
      plotRes(scResults$rr.summary.all, "Population growth rate",
        obs = scResults$obs.all,
        lowBound = 0, simRange = scResults$sim.all
      )
    })

    # lambda
    output$plot5 <- renderPlot({
      scResults <- dataInput1()
      plotRes(scResults$rr.summary.all, "Female population size",
        obs = scResults$obs.all,
        lowBound = 0
      )
    })

    output$plot6 <- renderPlot({
      out <- dataInput()
      dist <- out$oo$simDisturbance
      dist %>%
        tidyr::pivot_longer(c(Anthro, fire_excl_anthro),
          names_to = "dist_type",
          values_to = "dist"
        ) %>%
        mutate(dist_type = ifelse(dist_type == "Anthro", "Anthropogenic",
          "Fire excluding anthropogenic"
        )) %>%
        ggplot2::ggplot(ggplot2::aes(Year, dist)) +
        ggplot2::geom_point() +
        ggplot2::facet_wrap(~dist_type, nrow = 1) +
        ggplot2::theme_classic() +
        ggplot2::xlab("Year") +
        ggplot2::ylab("% disturbance") +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(size = 14),
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 14),
          axis.title.x = ggplot2::element_text(size = 16, face = "bold"),
          axis.title.y = ggplot2::element_text(size = 16, face = "bold"),
          strip.text = ggplot2::element_text(size = 14),
          strip.background = ggplot2::element_blank()
        ) +
        ggplot2::ylim(0, 100)
    })

    output$plot7 <- renderPlot({
      scResults <- dataInput1()
      plotRes(scResults$ksDists, "Adult female survival")
    })
    output$plot8 <- renderPlot({
      scResults <- dataInput1()
      plotRes(scResults$ksDists, "Recruitment")
    })
    output$plot9 <- renderPlot({
      scResults <- dataInput1()
      plotRes(scResults$ksDists, "Population growth rate")
    })

    # Download output #######
    observeEvent(input$saveOutput, {
      if (!is.integer(input$saveOutput)) {
        path <- parseSavePath(volumes(), input$saveOutput)$datapath
        if (grepl("\\.csv$", path)) {
          write.csv(dataInput1()$rr.summary.all, path, row.names = FALSE)
        } else if (grepl("\\.rds$", path)) {
          saveRDS(dataInput()$result, path)
        } else {
          warning(
            "file ", path, " could not be saved. Ensure that the file",
            " name includes .csv or .rds as the extension."
          )
        }
      }
    })

    # report.rmd is missing
    # output$downloadOutput <- downloadHandler(
    #   filename = function() {
    #     paste(input$inputIdCSV,".",sep = '', switch(
    #       input$format, HTML = 'html', PDF = 'pdf',  Word = 'docx'))
    #
    #
    #   },
    #
    #   content = function(file) {
    #     src <- normalizePath('report.Rmd')
    #
    #     # temporarily switch to the temp dir, in case you do not have write
    #     # permission to the current working directory
    #     owd <- setwd(tempdir())
    #     on.exit(setwd(owd))
    #     # on.exit(setwd(wdir)) # try this if download button won't work (write to working directory)
    #     file.copy(src, 'report.Rmd', overwrite = TRUE)
    #
    #     out <- render('report.Rmd', switch(
    #       input$format,
    #       PDF = pdf_document(), HTML = html_document(), Word = word_document()
    #     ))
    #     file.rename(out, file)
    #   }
    # )
  }

  # Run the application ---------------------------
  shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
}
