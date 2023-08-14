#' button.start
panel_about <- function(){
  div(class="about",
      div(HTML("&nbsp;")),
      div(HTML("&nbsp;")),
      fluidRow(
        column(12,
               column(12,
                      tags$h3("About GrowthTracker"),
                      HTML("<p>This app processes data from a CSV or an Excel file ...<br>
                      Tecnical details about the methods used here can be found in <b><a href='https://www.medrxiv.org/content/10.1101/2022.01.01.21268131v2'>this research article</a></b>.
                    </p>
                    "),
                      actionButton("button.start", "Start here"),
               ),
               column(12,
                      tags$h3("Got any feedback?"),
                      HTML("<p>...</p>")
               )
        )
      )
  )
}

panel_start_title <- function(){
  titlePanel(title = div(h3(HTML('Welcome to<br/>GrowthTracker'),
                            style = "padding-left: 15px; color: #000000; font-size: 25px; weigth:bold"),
                         h4(HTML(paste("Upload and explore your data, estimate",
                                       "growth rates, create short-term forecast",
                                       "and save your work with GrowthTracker",
                                       sep = "<br/>")),
                            style = "padding-left: 15px; padding-bottom: 15px; color: #797979; font-size: 15px; weigth:bold")),
             windowTitle = "Welcome")
}

#' start.button.new start.button.load
panel_start_main <- function(){
  fillPage(
    fluidRow(),
    #div(style="display:inline-block",submitButton("Analysis"))
    #div(style="display:inline-block",downloadButton('downloadData', 'Download Data'))
    #fluidRow(width = 12,
    #         style = "border-style: solid; border-color: black",
    #         align = "center",
    #         column(12, align = "center", "How to center this?")
    #         ),
    fluidRow(
      column(12, align = "left", #"center",
             style = 'padding-left:50px; padding-right:1px; padding-top:50px; padding-bottom:5px',
             #div(style="display:inline-block", submitButton("start.button.new", text = "Start new project"),
             actionButton(inputId = "start.button.new",
                          label = "Start new project",
                          style = "height: 60px; width: 180px; background-color:#9FB9B5; border-color:#9FB9B5")
             )
      ),
    fluidRow(
      column(12, align = "left", #"center",
             style = 'padding-left:50px; padding-right:1px; padding-top:5px; padding-bottom:10px',
             #div(style="display:inline-block", submitButton("start.button.load", text = "Load project"))
             actionButton(inputId = "start.button.load",
                          label = "Load existing project",
                          style = "height: 60px; width: 180px; background-color:#9FB9B5; border-color:#9FB9B5")
             )
      )
  )
}

panel_dashboard_title <- function(){
  titlePanel(title = div(h3(HTML('GrowthTracker dashboard'),
                            style = "padding-left: 15px; color: #000000; font-size: 25px; weigth:bold"),
                         h4(HTML(paste("Welcome to your dashboard.",
                                       "From here you can explore your data, run",
                                       "models and review nowcasting and",
                                       "forecasting based on your data",
                                       sep = "<br/>")),
                            style = "padding-left: 15px; padding-bottom: 15px; color: #797979; font-size: 15px; weigth:bold")),
             windowTitle = "Dashboard")
}

# dash.button.load dash.button.data dash.button.growth dash.button.nowcasting dash.button.forecasting
# dash.text.name
panel_dashboard_main <- function(){
  fillPage(
    fluidRow(),
    fluidRow(
      column(1),
      column(1,
             style = 'padding-left:1px; padding-right:1px; padding-top:50px; padding-bottom:5px',
             actionButton(inputId = "dash.button.load",
                          label = HTML("Load your<br/>dataset"),
                          style = "height: 120px; width: 120px; background-color:#9FB9B5; border-color:#9FB9B5")),
      column(1),
      column(1,
             style = 'padding-left:1px; padding-right:1px; padding-top:50px; padding-bottom:5px',
             actionButton(inputId = "dash.button.data",
                          label = HTML("Review your<br/>data"),
                          style = "height: 120px; width: 120px; background-color:#9FB9B5; border-color:#9FB9B5")),
      column(1),
      column(1,
             style = 'padding-left:1px; padding-right:1px; padding-top:50px; padding-bottom:5px',
             actionButton(inputId = "dash.button.growth",
                          label = HTML("Calculate<br/>growth rate<br/>estimates"),
                          style = "height: 120px; width: 120px; background-color:#9FB9B5; border-color:#9FB9B5")),
      column(1),
      column(1,
             style = 'padding-left:1px; padding-right:1px; padding-top:50px; padding-bottom:5px',
             actionButton(inputId = "dash.button.nowcasting",
                          label = HTML("Review<br/>nowcasting"),
                          style = "height: 120px; width: 120px; background-color:#9FB9B5; border-color:#9FB9B5")),
      column(1),
      column(1,
             style = 'padding-left:1px; padding-right:1px; padding-top:50px; padding-bottom:5px',
             actionButton(inputId = "dash.button.forecasting",
                          label = HTML("Review<br/>forecasting"),
                          style = "height: 120px; width: 120px; background-color:#9FB9B5; border-color:#9FB9B5")),
      column(1)
    ),
    fluidRow(column(1),
             column(11,
                    tags$div(id = "inline",
                             textInput(inputId = "dash.text.name",
                                       label = HTML("Project name: &nbsp;"),
                                       value = "New project")
                             )
                    ),
             style = "padding-top: 30px")
  )
}

#' 
panel_load_heading <- function(){
  titlePanel(title = div(h3('Load your dataset (step 1)',
                            style = "padding-left: 15px; color: #000000; font-size: 25px; weigth:bold"),
                         h4(HTML(paste("Load your file. If the data displayed does",
                                       "not coincide with your data, change the",
                                       "settings and try again",
                                       sep = "<br/>")),
                            style = "padding-left: 15px; padding-bottom: 15px; color: #797979; font-size: 15px; weigth:bold")),
             #title = "Load your dataset (step 1)",
             windowTitle = "Load your dataset")
}

#' load.browse.file
#' load.text.accept load.text.settings load.text.cancel
#' load.button.accept load.button.settings load.button.cancel
panel_load_sidebar <- function(){
  sidebarPanel(
    fileInput(
      #width = 4,
      inputId = "load.browse.file", # csvFile
      label = "Drag and drop the .csv or .xlsx file you want to process.\n",
      multiple = FALSE,
      buttonLabel = "Browse...",
      placeholder = "No file selected",
      accept = c(".csv", ".xls", "xlsx"),
    ),
    tags$span(
      #style = "vertical-align: bottom;",
      #actionButton("load.question", "", icon = icon("question")),
      p(id = "load.text.accept","If the data is correct, press", strong("Next"), style = "font-size:12px;"),
      p(id = "load.text.settings","If the data is not correct, press", strong("Settings"), style = "font-size:12px;"),
      p(id = "load.text.cancel","If you want to go to the previous page, press", strong("Cancel"), style = "font-size:12px;"),
      actionButton("load.button.accept", "Next"), # Upload data to project
      actionButton("load.button.settings", "Settings"),
      actionButton("load.button.cancel", "Cancel")
    )
  )
}

# load.table.data
panel_load_main <- function(){
  #mainPanel(fluidRow(dataTableOutput('load.table.data')))
  mainPanel(
    fluidRow(
      column(
        width = 12,
        dataTableOutput(outputId = "load.table.data")
        )
    )
  )
}

#' 
panel_load2_heading <- function(){
  titlePanel(title = div(h3('Check your dataset (step 2)',
                            style = "padding-left: 15px; color: #000000; font-size: 25px; weigth:bold"),
                         h4(HTML(paste("The file must have these columns:",
                                       "- Date",
                                       "- Count of positives",
                                       "- Count of tests",
                                       "and one optional column:",
                                       "- Category",
                                       sep = "<br/>")),
                            style = "padding-left: 15px; padding-bottom: 15px; color: #797979; font-size: 15px; weigth:bold")),
             #title = "Load your dataset (step 1)",
             windowTitle = "Load your dataset")
}

# load2.text.accept load2.text.cancel
panel_load2_sidebar <- function(){
  sidebarPanel(#tags$span(
    #style = "vertical-align: bottom;",
    #actionButton("load.question", "", icon = icon("question")),
    p(id = "load2.text.accept","When finished, press", strong("Next"), "(all columns must be different)", style = "font-size:12px;"),
    p(id = "load2.text.cancel","If you want to go to the previous page, press", strong("Cancel"), style = "font-size:12px;"),
    actionButton("load2.button.accept", "Next"),
    actionButton("load2.button.cancel", "Cancel")
  )
}

# load2.select.date load2.select.dateformat load2.select.positives load2.select.tests load2.select.categories
# load2.text.date
panel_load2_main <- function(){
  #print(colnames(dataInput$data))
  column(12,
         style='height: 400px; overflow-y: scroll;',#style='scroll; height: 400px; overflow-y: scroll;',
       fluidRow(),
      fluidRow(column(width = 1),
               column(width = 10, 
                      style = "padding-right: 15px; border: 1px solid; border-color: #9FB9B5",
                      fluidRow(
                        column(width = 4,
                               selectInput(inputId = "load2.select.date",
                                           label = div("Which column has the date/time?",
                                                       style = "padding-bottom: 5px; padding-top: 5px;"),
                                           #label = div(p("Which column has the date?",
                                          #                style = "padding-left: 5px; padding-bottom: 5px; color: #000000; font-size: 15px; weigth:normal")),
                                           choices = list(defaultColumn),
                                           selected = defaultColumn),
                               selectInput(inputId = "load2.select.dateformat",
                                         label = "Indicate date format:",
                                         choices = list(defaultSelectFormat),
                                         selected = defaultSelectFormat),
                               style = "border-right: 1px dashed #9FB9B5"
                        ),#padding-left: 15px;
                        column(width = 8,
                               div(tableOutput(outputId = "load2.text.date"), style = "font-size:80%; padding-top: 15px"),
                               #wellPanel(div(tableOutput(outputId = "load2.text.date"), style = "font-size:80%"),
                                                    #htmlOutput(outputId = "load2.text.date"),
                                #                    style = "padding-top: 15px; padding-left: 15px; color: #000000; font-size: 15px")
                                                    #div(h5('Hola',# TODO
                                                    #       style = "padding-left: 15px; color: #000000; font-size: 15px")))
                        )
                      ),
               ),
               column(width = 1),
               style = "padding-bottom: 2px"
      ),
      fluidRow(column(width = 1),
               column(width = 10, 
                      style = "padding-right: 15px; border: 1px solid; border-color: #9FB9B5",
                      fluidRow(
                        column(width = 4,
                               selectInput(inputId = "load2.select.positives",
                                           label = div("Which column has the positives?",
                                                       style = "padding-bottom: 5px; padding-top: 5px;"),
                                           #label = div(h5("Which column has the positives?",
                                           # style = "padding-left: 5px; padding-bottom: 5px; color: #000000; font-size: 15px; weigth:normal")),
                                           choices = list(defaultColumn),
                                           selected = defaultColumn),
                               style = "border-right: 1px dashed #9FB9B5"
                        ),#padding-left: 15px; 
                        column(width = 8, div(h5('',
                                                 style = "padding-left: 15px; color: #000000; font-size: 15px"))
                        )
                      ),
               ),
               column(width = 1),
               style = "padding-bottom: 2px"
      ),
      fluidRow(column(width = 1),
               column(width = 10, 
                      style = "padding-right: 15px; border: 1px solid; border-color: #9FB9B5",
                      fluidRow(
                        column(width = 4, selectInput(inputId = "load2.select.tests",
                                                      label = div("Which column has the number of tests?",
                                                                  style = "padding-bottom: 5px; padding-top: 5px;"),
                                                      #label = div(h5("Which column has the number of tests?",
                                                      #               style = "padding-left: 5px; padding-bottom: 5px; color: #000000; font-size: 15px; weigth:normal")),
                                                      choices = list(defaultColumn),
                                                      selected = defaultColumn),
                               style = "border-right: 1px dashed #9FB9B5"
                        ),#padding-left: 15px; 
                        column(width = 8, div(h5('',
                                                 style = "padding-left: 15px; color: #000000; font-size: 15px"))
                        )
                      ),
               ),
               column(width = 1),
               style = "padding-bottom: 2px"
      ),
      fluidRow(column(width = 1),
               column(width = 10, 
                      style = "padding-right: 15px; border: 1px solid; border-color: #9FB9B5",
                      fluidRow(
                        column(width = 4, selectInput(inputId = "load2.select.categories",
                                                      label = div("(Optional) Which column has the list of categories?",
                                                                  style = "padding-bottom: 5px; padding-top: 5px;"),
                                                      #label = div(h5("(Optional) Which column has the list of categories?",
                                                      #               style = "padding-left: 5px; padding-bottom: 5px; color: #000000; font-size: 15px; weigth:normal")),
                                                      choices = list(defaultColumn),
                                                      selected = defaultColumn),
                               style = "border-right: 1px dashed #9FB9B5"
                        ),#padding-left: 15px; 
                        column(width = 8, div(h5('',
                                                 style = "padding-left: 15px; color: #000000; font-size: 15px"))
                        )
                      ),
               ),
               column(width = 1)
      )
    )
}

# ----
debug_message <- function(title = "", message = "", error = ""){
  showModal(modalDialog(
    title = title,
    span(message),
    div(span("Details: ", code(error), ".")),
    easyClose = TRUE,
    footer = modalButton("OK")
  ))
}


