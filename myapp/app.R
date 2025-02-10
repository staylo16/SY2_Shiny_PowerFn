#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#Published: https://uoe-maths.shinyapps.io/SY2-EstimatorTheory/





z_test_power <- function(mu_star_case = 0, 
                         alt = c("neq","<",">"),
                         mu_0 = 0, sigma2 = 1, n = 10,
                         alpha = 0.05, plot = TRUE, show_all = FALSE){
 
 mu_star <- mu_0 + seq(-3.5, 3.5, by=0.01)
 
 mu_star_case <- max(mu_star_case, min(mu_star))
 mu_star_case <- min(mu_star_case, max(mu_star))
 alt <- match.arg(alt,c("neq","<",">"))
 mu_0 <- max(mu_0, min(mu_star))
 mu_0 <- min(mu_0, max(mu_star))
 sigma2 <- max(sigma2, sqrt(.Machine$double.eps))
 n <- ceiling(max(n,2))
 alpha <- max(alpha,0)
 alpha <- min(alpha,1)
 
 
 if(alt=="neq"){
  POWER_VALUE <- pnorm(-qnorm(1-alpha/2)-(mu_star_case - mu_0)/sqrt(sigma2/n), mean = 0, sd = 1) + 
   1 - pnorm(qnorm(1-alpha/2)-(mu_star_case - mu_0)/sqrt(sigma2/n), mean = 0, sd = 1)
 }else if(alt=="<"){
  POWER_VALUE <- pnorm(-qnorm(1-alpha)-(mu_star_case - mu_0)/sqrt(sigma2/n), mean = 0, sd = 1)
 }else{
  POWER_VALUE <- 1 - pnorm(qnorm(1-alpha)-(mu_star_case - mu_0)/sqrt(sigma2/n), mean = 0, sd = 1)
 }
 if(!plot){return(POWER_VALUE)}
 
 #mu_0 <- 0
 #mu_star_case <- 0.5
 #alt <- c("neq","<",">")[3]
 #sigma2 <- 4
 #n <- 30
 #alpha <- 0.1
 
 pow_2 <- pnorm(-qnorm(1-alpha/2)-(mu_star - mu_0)/sqrt(sigma2/n), mean = 0, sd = 1) + 
  1 - pnorm(qnorm(1-alpha/2)-(mu_star - mu_0)/sqrt(sigma2/n), mean = 0, sd = 1) 
 pow_low <- pnorm(-qnorm(1-alpha)-(mu_star - mu_0)/sqrt(sigma2/n), mean = 0, sd = 1) 
 pow_upp <- 1 - pnorm(qnorm(1-alpha)-(mu_star - mu_0)/sqrt(sigma2/n), mean = 0, sd = 1) 
 
 par(mfrow=c(2,1))
 plot(mu_star, pow_2, type="n",ylim=c(0,1), ylab="Power", main = "Power Curve, Z-test",
      xlab = expression(mu^"*"),cex.lab = 1.2)
 grid()
 abline(h = alpha, lty=2)
 if(show_all | alt == "neq") lines(mu_star, pow_2, col=2,lwd=2)
 if(show_all | alt == "<") lines(mu_star, pow_low, col=3,lwd=2)
 if(show_all | alt == ">") lines(mu_star, pow_upp, col=4,lwd=2)
 points(mu_star_case,POWER_VALUE,pch=16,cex=2) 
 
 #mu_star_case <- sample(mu_star,size=1)
 
 x <- seq(-8,8,by=0.01)
 z_h0 <- dnorm(x,mean=0,sd=1)
 z_h1 <- dnorm(x - (mu_star_case-mu_0)/sqrt(sigma2/n),mean=0,sd=1)
 
 plot(NA,NA,xlab="z",ylab="f(z)",xlim=range(x),ylim=c(0,max(z_h0,z_h1)), 
      main = "Distribution of the Test Statistic")
 if(alt=="neq"){
  z1 <- qnorm(alpha/2)
  z2 <- qnorm(1-alpha/2)
  polygon(c(min(x),x[x<z1],z1,z1),
          c(0,z_h0[x<z1],dnorm(z1),0), col=2,density = 20,angle = -45)
  polygon(c(z2,z2,x[x>z2],max(x)),
          c(0,dnorm(z2),z_h0[x>z2],0),col=2,density = 20,angle = -45)
  polygon(c(min(x),x[x<z1],z1,z1),
          c(0,z_h1[x<z1],dnorm(z1 - (mu_star_case-mu_0)/sqrt(sigma2/n)),0), col=3,density=20)
  polygon(c(z2,z2,x[x>z2],max(x)),
          c(0,dnorm(z2 - (mu_star_case-mu_0)/sqrt(sigma2/n)),z_h1[x>z2],0),col=3,density = 20)
  abline(v=c(z1,z2),lty=2,lwd=2,col=1)
  lines(x,z_h1,lwd=2,col=2)
  
 }else if(alt=="<"){
  z<- qnorm(alpha)
  polygon(c(min(x),x[x<z],z,z),
          c(0,z_h0[x<z],dnorm(z),0), col=2,density = 20,angle = -45)
  polygon(c(min(x),x[x<z],z,z),
          c(0,z_h1[x<z],dnorm(z - (mu_star_case-mu_0)/sqrt(sigma2/n)),0), col=3,density=20)
  abline(v=z,lty=2,lwd=2,col=1)
  lines(x,z_h1,lwd=2,col=3)
  
 }else{
  z <- qnorm(1-alpha)
  polygon(c(z,z,x[x>z],max(x)),
          c(0,dnorm(z),z_h0[x>z],0),col=2,density = 20,angle = -45)
  polygon(c(z,z,x[x>z],max(x)),
          c(0,dnorm(z-(mu_star_case-mu_0)/sqrt(sigma2/n)),z_h1[x>z],0),col=3,density = 20)
  abline(v=z,lty=2,lwd=2,col=1)
  lines(x,z_h1,lwd=2,col=4)
 }
 
 lines(x,z_h0,lwd=2)
 abline(v = 0,lty=2)
 abline(v = (mu_star_case-mu_0)/sqrt(sigma2/n),lty=2)
 
 
 return(invisible(POWER_VALUE))
 
}



one_sample_t_test_power <- function(mu_star_case = 0, 
                                    alt = c("neq","<",">"),
                                    mu_0 = 0, s2 = 1, n = 10,
                                    alpha = 0.05, plot = TRUE, show_all = FALSE){
 
 mu_star <- mu_0 + seq(-3.5, 3.5, by=0.01)
 
 
 mu_star_case <- max(mu_star_case, min(mu_star))
 mu_star_case <- min(mu_star_case, max(mu_star))
 alt <- match.arg(alt,c("neq","<",">"))
 mu_0 <- max(mu_0, min(mu_star))
 mu_0 <- min(mu_0, max(mu_star))
 s2 <- max(s2, sqrt(.Machine$double.eps))
 n <- ceiling(max(n,2))
 alpha <- max(alpha,0)
 alpha <- min(alpha,1)
 
 
 if(alt=="neq"){
  POWER_VALUE <- pt(-qt(1-alpha/2, df=n-1)-(mu_star_case - mu_0)/sqrt(s2/n), df=n-1) + 
   1 - pt(qt(1-alpha/2, df=n-1)-(mu_star_case - mu_0)/sqrt(s2/n), df=n-1)
 }else if(alt=="<"){
  POWER_VALUE <- pt(-qt(1-alpha, df=n-1)-(mu_star_case - mu_0)/sqrt(s2/n), df=n-1)
 }else{
  POWER_VALUE <- 1 - pt(qt(1-alpha, df=n-1)-(mu_star_case - mu_0)/sqrt(s2/n), df=n-1)
 }
 if(!plot) return(POWER_VALUE)
 
 
 pow_2 <- pt(-qt(1-alpha/2, df=n-1)-(mu_star - mu_0)/sqrt(s2/n), df=n-1) + 
  1 - pt(qt(1-alpha/2, df=n-1)-(mu_star - mu_0)/sqrt(s2/n), df=n-1) 
 pow_low <- pt(-qt(1-alpha, df=n-1)-(mu_star - mu_0)/sqrt(s2/n), df=n-1) 
 pow_upp <- 1 - pt(qt(1-alpha, df=n-1)-(mu_star - mu_0)/sqrt(s2/n), df=n-1) 
 
 par(mfrow=c(2,1))
 plot(mu_star, pow_2, type="n",ylim=c(0,1), ylab="Power", main = "Power Curve, one-sample t-test",
      xlab = expression(mu^"*"),cex.lab = 1.2) 
 grid()
 abline(h = alpha, lty=2)
 if(show_all | alt == "neq")lines(mu_star, pow_2, col=2,lwd=2)
 if(show_all | alt == "<")lines(mu_star, pow_low, col=3,lwd=2)
 if(show_all | alt == ">")lines(mu_star, pow_upp, col=4,lwd=2)
 points(mu_star_case,POWER_VALUE,pch=16,cex=2) 
 
 
 x <- seq(-8,8,by=0.01)
 z_h0 <- dt(x,df = n-1)
 z_h1 <- dt(x-(mu_star_case-mu_0)/sqrt(s2/n), df=n-1)
 
 plot(NA,NA,xlab="t",ylab="f(t)",xlim=range(x),ylim=c(0,max(z_h0,z_h1)), 
      main = "Distribution of the Test Statistic")
 if(alt=="neq"){
  z1 <- qt(alpha/2,n-1)
  z2 <- qt(1-alpha/2,df=n-1)
  polygon(c(min(x),x[x<z1],z1,z1),
          c(0,z_h0[x<z1],dt(z1,df=n-1),0), col=2,density = 20,angle = -45)
  polygon(c(z2,z2,x[x>z2],max(x)),
          c(0,dt(z2,df=n-1),z_h0[x>z2],0),col=2,density = 20,angle = -45)
  polygon(c(min(x),x[x<z1],z1,z1),
          c(0,z_h1[x<z1],dt(z1-(mu_star_case-mu_0)/sqrt(s2/n),df=n-1),0), col=3,density=20)
  polygon(c(z2,z2,x[x>z2],max(x)),
          c(0,dt(z2-(mu_star_case-mu_0)/sqrt(s2/n),df=n-1),z_h1[x>z2],0),col=3,density = 20)
  abline(v=c(z1,z2),lty=2,lwd=2,col=1)
  lines(x,z_h1,lwd=2,col=2)
  
 }else if(alt=="<"){
  z<- qt(alpha,df=n-1)
  polygon(c(min(x),x[x<z],z,z),
          c(0,z_h0[x<z],dt(z,df=n-1),0), col=2,density = 20,angle = -45)
  polygon(c(min(x),x[x<z],z,z),
          c(0,z_h1[x<z],dt(z-(mu_star_case-mu_0)/sqrt(s2/n),df=n-1),0), col=3,density=20)
  abline(v=z,lty=2,lwd=2,col=1)
  lines(x,z_h1,lwd=2,col=3)
  
 }else{
  z <- qt(1-alpha,df=n-1)
  polygon(c(z,z,x[x>z],max(x)),
          c(0,dt(z,df=n-1),z_h0[x>z],0),col=2,density = 20,angle = -45)
  polygon(c(z,z,x[x>z],max(x)),
          c(0,dt(z-(mu_star_case-mu_0)/sqrt(s2/n),df=n-1),z_h1[x>z],0),col=3,density = 20)
  abline(v=z,lty=2,lwd=2,col=1)
  lines(x,z_h1,lwd=2,col=4)
 }
 
 lines(x,z_h0,lwd=2)
 abline(v = mu_0,lty=2)
 abline(v = (mu_star_case-mu_0)/sqrt(s2/n),lty=2)
 
 return(invisible(POWER_VALUE))
}




library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
 
 
 # Application title
 titlePanel(withMathJax("Power Curves")),
 
 # Sidebar with a slider input for number of bins 
 sidebarLayout(
  sidebarPanel(
   
   
   selectInput("test", "Hypothesis Test:", c("Z-test" = "z",
                                             "One-sample t-test" = "t")),
   selectInput("alt", "Alternative:", c("two-sided" = "neq",
                                        "one-sided, less than" = "<",
                                        "one-sided, more than" = ">")),
   sliderInput("alpha", "Significance Level", min=0.01, max=0.5, value = 0.05,step =0.01),
   sliderInput("n", "Sample Size", min=3, max=100, value=10, step= 1),
   sliderInput("sig2", "Population/Sample Variance", min=1, max=16, value=1, step=0.5),
   sliderInput("mustar", "Mu-star, case under H1", min=-3, max=3, value=0.9, step=0.1),
   
   
   h3("Power:"),
   textOutput("Power",inline=TRUE),
   
   h3("Key:"),
   h4("{YET TO BE COMPLETED!!!}"),
   
   width = 3
   
  ),
  
  
  # Show a plot of the generated distribution
  mainPanel(
   plotOutput("Plots", height = "1000px")
  )
 )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
 r <- reactiveValues(seed = as.numeric(Sys.time()))
 
 
 output$Plots <- renderPlot({
  
  if(input$test == "z"){
   pow_val <- z_test_power(mu_star_case = as.numeric(input$mustar), 
                           alt = input$alt, 
                           mu_0 = 0, 
                           sigma2 = as.numeric(input$sig2), 
                           n = as.numeric(input$n), 
                           alpha = as.numeric(input$alpha))
  }else{
   pow_val <- one_sample_t_test_power(mu_star_case = as.numeric(input$mustar), 
                                      alt = input$alt, 
                                      mu_0 = 0, 
                                      s2 = as.numeric(input$sig2), 
                                      n = as.numeric(input$n), 
                                      alpha = as.numeric(input$alpha))
  }
  output$Power <- renderText({round(pow_val, 4)})
 })
 
 
 
}

# Run the application 
shinyApp(ui = ui, server = server)
