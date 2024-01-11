library(shiny)
library(data.table)
library(shinyTime)
library(markdown)

transfer_csv.num.col <- function(csv.file){
  numcol <- suppressWarnings(as.numeric(gsub("X", "",  colnames(csv.file))))
  wavelength <-  numcol[suppressWarnings(which( !is.na(numcol) & numcol > 100))]
  numcol <- suppressWarnings(which( !is.na(numcol) & numcol > 100))
  
  returnlist <- list(numcol, wavelength)
  names(returnlist) <- c("numcol", "wl")
  return(returnlist)
}

xaxisdate <- function (dayvec, type = NA, formatd = "day", ydata = NA, pch = 20, 
                       cex = 1, coltypep = "red", las = 2, tz = "UTC") 
{
  dayvec <- as.POSIXct(dayvec, tz = tz, origin = "1970-01-01")
  if (is.na(type)) 
    type <- "n"
  dayvecc <- as.Date(dayvec, tz = tz)
  datep <- c()
  datepp <- unique(dayvecc)
  for (i in 1:length(datepp)) {
    datep[i] <- which(dayvecc == datepp[i])[1]
  }
  if (formatd == "cw") {
    datepp <- lubridate::week(datepp)
    datep <- datep[which(!duplicated(datepp))]
    datepp <- datepp[which(!duplicated(datepp))]
  }
  if (formatd == "time") {
    datepp <- dayvec
    datep <- which(!duplicated(lubridate::hour(dayvec)))
    datepp <- datepp[which(!duplicated(lubridate::hour(dayvec)))]
    datepp <- substr(datepp, 12, 16)
  }
  if (formatd == "day") {
    datepp <- datepp
    datepp <- paste(lubridate::month(datepp, label = T, abbr = T), 
                    lubridate::day(datepp), sep = "-")
  }
  if (formatd == "year-month") {
    datepp <- strftime(datepp, format = "%y-%m", tz = tz)
    datep <- datep[which(!duplicated(datepp))]
    datepp <- datepp[which(!duplicated(datepp))]
  }
  padj <- 0.5
  axis(1, at = datep, datepp, padj = padj, las = las, cex.axis = cex)
  axis(2)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       lwd = 1.2)
  if (type == "l") 
    abline(v = datep, lty = 3)
  if (type == "p") 
    points(datep, ydata[datep], pch = pch, col = coltypep, 
           cex = cex)
  box()
}

options(shiny.maxRequestSize=1500*1024^2)

# UI Definition
ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      /* Adjust the font size of input labels */
      .control-label {
        font-size: 12px;  /* Adjust this value as needed */
      }
      /* Adjust the font size within input boxes */
      input, select, textarea {
        font-size: 12px;  /* Adjust this value as needed */
        height: auto;     /* Adjust height as needed */
        padding: 5px;     /* Adjust padding for better text visibility */
      }
    "))
  ),
  
  titlePanel(title=div(img(src="dt_logo.png", height="5%", width="5%"), "DT:SVD Analysis of Spectral Data")
             , windowTitle = "DT:SVD Analysis of Spectral Data"),
  
  # Adjusting the sidebar layout to be narrower
  sidebarLayout(
    sidebarPanel(
      width = 3,  # Adjust the width of the sidebar
      
      fileInput("file_upload", "", multiple = T, accept = c(".csv"), buttonLabel = "csv-Dateien", placeholder = NA, width = "100%"),
      checkboxInput("meancenter", "Mean Center Data", TRUE),
      numericInput("rank", "Choose Rank", 1, min = 1, max = 409),
      
      sliderInput("xrange.rank", "Choose Rank Range (zoom)", min = 1, max = 409, value = c(1, 409)),
      sliderInput("xrange.lambda", "Choose lambda for SVD (new calculation)", min = 190, max = 598, value = c(200, 588)),
      
      uiOutput("dateTimeRangeInput")  # Correctly placed within sidebarPanel
    ),
    
    mainPanel(
      tabsetPanel(
        type = "tabs",
        
        tabPanel("SVD Analysis", 
                 plotOutput("plot1"),
                 fluidRow(  # Placing plot2 and plot3 next to each other
                   column(width = 6, plotOutput("plot2")),
                   column(width = 6, plotOutput("plot3"))
                 )
        ),
        
        tabPanel("Über",
                 fluidRow(
                   column(6,
                          includeMarkdown("about.Rmd")# Markdown("about.md")
                   )
                 )
        ),
        
        tabPanel("Versionsverlauf",
                 fluidRow(
                   column(6,
                          includeMarkdown("Version.Rmd")# Markdown("about.md")
                   )
                 )
        ),
        
        tabPanel("Github",
                 fluidRow(
                   column(6,
                          includeMarkdown("Github.Rmd")# Markdown("about.md")
                   )
                 )
        )
      )
    )
  )
)

# Server Logic
server <- function(input, output) {
  
  # Reactive expression to read and process the data
  data_reactive <- reactive({
    req(input$file_upload)
    
    # Read the file
    df <- fread(input$file_upload$datapath, sep = ";", dec = ",")
    df <- df[ order(df$datetime) , ]
    df$datetime <- as.POSIXct( as.character( df$datetime ), tz = "UTC" )
    
    # Perform operations
    ppp <- transfer_csv.num.col(df)
    lambda.start <- ppp$wl
    # subset df by input$xrange.lambda
    xrange.lambda.sub <- which(ppp$wl %in% input$xrange.lambda)
    
    xrange.df.sub <- c(which( !(1 : ncol( df )) %in% ppp$numcol)
                       , ppp$numcol[ xrange.lambda.sub[ 1 ] : xrange.lambda.sub[ 2 ]])
    
    df <- df[ , xrange.df.sub, with = F]
    ppp <- transfer_csv.num.col(df)
    
    datetime <- df$datetime
    lambda <- ppp$wl
    spc <- df[, ppp$numcol, with = FALSE]
    
    # Mean centering if selected
    if (input$meancenter) {
      mean_spc <- spc[, lapply(.SD, mean, na.rm = TRUE)]
      spc <- t(apply(spc, 1, function(x) x - as.numeric(mean_spc)))
    }
    
    svd_result <- svd(spc)
    list(datetime = datetime, lambda = lambda, svd = svd_result, lambda.start = lambda.start)
  })
  
  # Update slider ####
  observe({
    data <- data_reactive()
    if (!is.null(data$lambda)) {
      updateSliderInput(getDefaultReactiveDomain(), "xrange.rank", min = 1, max = length(data$lambda.start))
      updateSliderInput(getDefaultReactiveDomain(), "xrange.lambda", min = min(data$lambda.start), max = max(data$lambda.start))
    }
  })
  
  # Datetime
  output$dateTimeRangeInput <- renderUI({
    data <- data_reactive()
    
    if (!is.null(data$datetime)) {
      minDate <- min(data$datetime, na.rm = TRUE)
      maxDate <- max(data$datetime, na.rm = TRUE)
      
      fluidRow(
        column(6, dateInput("startDate", "Start Date", value = minDate, min = minDate, max = maxDate)),
        column(6, timeInput("startTime", "Start Time", value = format(minDate, "%H:%M:%S"))),
        column(6, dateInput("endDate", "End Date", value = maxDate, min = minDate, max = maxDate)),
        column(6, timeInput("endTime", "End Time", value = format(maxDate, "%H:%M:%S"))),
        column(6, checkboxInput("xdatetime", "time-resolved or consecutively?", TRUE))
      )
    }
  })
  
  # # Generate datetime slider UI based on the uploaded file
  # output$datetimeSlider <- renderUI({
  #   data <- data_reactive()
  #   if (!is.null(data$datetime) && length(data$datetime) > 0 && !all(is.na(data$datetime))) {
  #     sliderInput("datetimeRange", "Select DateTime Range",
  #                 min = min(data$datetime, na.rm = TRUE),
  #                 max = max(data$datetime, na.rm = TRUE),
  #                 value = c(min(data$datetime, na.rm = TRUE), max(data$datetime, na.rm = TRUE)),
  #                 timeFormat = "%Y-%m-%d %H:%M:%S")
  #   }
  # })
  
  # Plot outputs
  output$plot1 <- renderPlot({
    
    data <- data_reactive()
    
    xrange.rank <- input$xrange.rank
    lengthp <- length(data$svd$d)
    xrangep.rank <- xrange.rank[ 1 ] : xrange.rank[ 2 ]
    
    colp <- rep("black", lengthp)
    colp[input$rank] <- "red"
    colp <- colp[xrangep.rank]
    
    cexp <- rep(1, lengthp)
    cexp[input$rank] <- 1.5
    cexp <- cexp[xrangep.rank]
    
    lwdp <- rep(1, lengthp)
    lwdp[input$rank] <- 2
    lwdp <- lwdp[xrangep.rank]
    
    x1 <- data$svd$d[xrangep.rank]
    
    par(mar = c(4, 4, 1, 1))
    plot(x1 / sqrt(nrow( data$svd$u)) 
         , type = 'p'
         , col = colp
         , cex = cexp
         , lwd = lwdp
         , log = "y"
         , xlab = "Rank", ylab = "Singular Value Decomposition"
         , main = "Singular Value Decomposition")
  })
  
  # Adjusting plot2 based on the datetime range selected
  output$plot2 <- renderPlot({
    data <- data_reactive()
    if (!is.null(data$datetime) && !is.null(input$startDate) && !is.null(input$endDate)) {
      startDate <- as.POSIXct(input$startDate)
      endDate <- as.POSIXct(input$endDate)
      startTime <- strftime(as.character(input$startTime), format = "%H:%M:%S")
      endTime <- strftime(as.character(input$endTime), format = "%H:%M:%S")
      
      filteredData <- data$datetime >= as.POSIXct( as.character( paste(startDate, startTime)), tz = "UTC") & 
        data$datetime <= as.POSIXct( as.character( paste(endDate, endTime)), tz = "UTC")
      
      if (any(filteredData)) {
        
        par(mar = c(4, 4, 3, 1))
        if(input$xdatetime) plot(data$datetime[filteredData], data$svd$u[filteredData, input$rank],
                                 xlab = "Date/Time", ylab = "Zeiteigenvektor"
                                 , main = paste("Zeiteigenvektor, Rank =", input$rank))
        
        if(!input$xdatetime){ plot(data$svd$u[filteredData, input$rank],
                                   xlab = "", ylab = "Zeiteigenvektor"
                                   , main = paste("Zeiteigenvektor, Rank =", input$rank)
                                   , axes = F)
          xaxisdate(data$datetime[filteredData])
        }
      }
    }
  })
  
  output$plot3 <- renderPlot({
    data <- data_reactive()
    
    xrange.rank <- input$xrange.rank
    xrangep.rank <- xrange.rank[ 1 ] : xrange.rank[ 2 ]
    lengthp <- length(data$svd$d)
    lambda <- data$lambda
    
    colp <- rep("black", lengthp)
    
    cexp <- rep(1, lengthp)
    
    lwdp <- rep(1, lengthp)
    
    x1 <- data$svd$v[, input$rank]
    
    par(mar = c(4, 4, 3, 1))
    plot(lambda, x1
         , type = 'p'
         , col = colp
         , cex = cexp
         , lwd = lwdp
         , xlab = expression(paste(lambda, " in nm")), ylab = "Spektrale Eigenvektoren"
         , main = "Spektrale Eigenvektoren")
  })
}

# Run the application
shinyApp(ui = ui, server = server)
