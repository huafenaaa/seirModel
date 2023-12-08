
library(shiny)
library ( shinydashboard)
library(deSolve)
library( cowplot)
library(ggplot2)
library (tidyverse)
library(ggrepel)
library (rstudioapi)
library ( showtext)
library(renv)
#设置工作路径
#current_dir=dirname(getSourceEditorContext( )$path)
#setwd( current_dir)
#设置中文显示支持
#showtext_auto()

#创建数学模型
seirModel <- function(time, state,parameters) {
  with(as.list(c(state, parameters)),{
    dSt = mu* ( 1 - S) - beta * S* I
    dEt= beta *S*I - mu *E - lamda*E
    dIt = - (mu+gamma)* I + lamda*E
    dRt= gamma *I -mu*R
    return(list(c( dSt,dIt,dRt,dEt)))
  })
}

#定义用户接口组件
ui <-dashboardPage(
  skin = "purple",
  dashboardHeader(title ="传染病模型仿真"),
  dashboardSidebar(
    tags$head(tags$style(HTML('
    .main-header  .logo {
     font-family: "Georgia",Times, "Times New Roman", serif;
     font-weight: bold;
     font-size:20px;
   }
'))),


    #定义用户滑块输入各参数信息
    width =300,
    sliderInput("mu",
                "出生率(死亡率) :",
                min=0, max=0.2,value =0.00
    ),
    sliderInput("ratioInfec",
                "初始感染比率:",
                min =0, max = 0.0002, ticks =FALSE,value=0.00002
    ),
    sliderInput("ratioSuscep",
                "初始易感比率:",
                min=0.8,max=1, ticks= FALSE,value= 0.9999
    ),
    sliderInput( "lamda",
                 "易感者转化为感染者的速率:",
                 min= 0,max =0.3, value=0.03
    ),
    sliderInput( "beta" ,
                 "疾病传播速率: ",
                 min=0,max=20,value= 4
    ),
    sliderInput( "recoveryInterval",
                 "感染者康复周期（天）: ",
                 min=1, max =47, value= 10
    ),
    sliderInput("timeinterval",
                "仿真时间〔周）:",
                min= 1,max= 156, value= 104
    )

  ),
  dashboardBody (
    #精绘制输出结果图形
    fluidRow(plotOutput ( "outcomePlot" )),
    br(),
    br(),
    br(),
    br(),
    fluidRow(
      #动态响应用户输入并动态输出
      valueBoxOutput("finalInfeciton", width= 6)
    ),
  )
)

#定义服务器端组件
server<- function(input,output) {
  #创建响应式输出
  reactiveData <- reactive({
    init    <-
      c(
        S=input$ratioSuscep,
        I = input$ratioInfec,
        R=0,
        E=1-input$ratioSuscep- input$ratioInfec
      )
    #仿真参数由用户界面输入
    parameters <-
      c(beta = input$beta,
        gamma =1/input$recoveryInterval,
        mu = input$mu,
        lamda = input$lamda)
    #评估时间
    times<-seq(0, input$timeinterval, by =1/7)
    ##求解一阶常微分方程
    outcome <- ode(
      y =init,
      times=times,
      func =seirModel,
      parms = parameters
    )
    #返问结果
    as.data.frame (outcome)
  })
  #渲染结果并绘制图形
  output$outcomePlot<- renderPlot({
    outcome<-
      reactiveData( )%>%
      pivot_longer(c(-time),
                   names_to = "key" , values_to = "value")%>%
      mutate(
        id= row_number(),
        感染状态= recode(
          key,
          S="易感者",
          I="感染者",
          R="康复者",
          E="潜伏者"
        ),
        leftSide = recode(
          key,
          S="易感者",
          I="感染者",
          R="",
          E=""
        ),
        rightSide= recode(
          key,
          S="",
          I="",
          R="康复者",
          E="潜伏者"
        )
      )

    #绘制图形
    ggplot(data = outcome,
           aes(
             x= time ,
             y= value,
             group=感染状态,
             col=感染状态,
             label =感染状态,
             data_id =id
           ))+
      ylab("人口比率")+xlab("时间(周)") +ylim(c(0,1))+
      geom_line(size = 2) +
      geom_text_repel(
        data = subset(outcome,time== max(time ) ) ,
        aes(label =rightSide),
        size = 7,
        segment.size =0.3,
        nudge_x = 0,
        hjust =0.5,vjust=0.5
      ) +
      geom_text_repel(
        data = subset(outcome, time == min(time) ) ,
        aes(label = leftSide),
        size = 7,
        segment.size=0.3,
        nudge_x = 0,
        hjust = 0.2,vjust=0.5
      )+
      scale_colour_manual(values=c("green ","red","purple","orange")) +
      scale_y_continuous( labels = scales::percent, limits = c(0,1))+theme_bw()

  })
  #统计最终感染人数占比
  output$finalInfeciton <- renderValueBox({
    valueBox("观察期内最终感染人群占比: ",
             reactiveData()%>% filter(time == max(time)) %>% select(R)%>%
               mutate(R =round(100* R,3))%>% paste0("%"),color="purple")
  })

}





# Run the application
shinyApp(ui = ui, server = server)
