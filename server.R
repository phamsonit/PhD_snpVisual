library(shiny)
library(shinythemes)
library(cluster)
library(shape)
library(ggplot2)
options(shiny.maxRequestSize = 50*1024^2)

function(input, output, session) {
  
  #load data from csv file
  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.

  #read patterns file
  read_pattern <- reactive({
    inFile <- input$datafile
    read_pattern <- read.csv(inFile$datapath,header = FALSE)
    #read_pattern <- read.csv("BC_pattern_sample.csv",header = FALSE)
  })
  
  #read index file
  read_index <- reactive({
    inFile <- input$indexfile
    read_index <- read.csv(inFile$datapath,header = FALSE)
    #read_index <- read.csv("rape/rapsodyn_filter_index.csv",header = FALSE)
    #read_index <- read.csv("BC_index.csv",header = FALSE)
    #read_index <- read.csv("chicken_SNP_index.csv",header = FALSE)
  })
  
  #read chromosome file
  read_chr <- reactive({
    inFile <- input$chrfile
    read_chr <- read.csv(inFile$datapath,header = FALSE)
    #read_chr <- read.csv("rape/rapsodyn_chr.csv",header = FALSE)
    #read_chr <- read.csv("BC_chr.csv",header = FALSE)
    #read_chr <- read.csv("chicken_chr.csv",header = FALSE)
  })
  
  #read gene file
  read_gene <- reactive({
    inFile <- input$genefile
    read_chr <- read.csv(inFile$datapath,header = FALSE)
    #read_gene <- read.csv("rape/rapsodyn_gene_info.csv",header = FALSE)
  })
  
  #read QTL file
  read_QTL <- reactive({
    inFile <- input$QTLfile
    read_chr <- read.csv(inFile$datapath,header = FALSE)
    #read_QTL <- read.csv("rape/rapsodyn_QTLs_Indiv.csv",header = FALSE)
    #read_QTL <- read.csv("chicken_QTL.csv",header = FALSE)
  })
  
  #create sliderInput to choose cutval
  output$cutval <- renderUI({
    if(!is.null(input$datafile))
    {
      m <- read_pattern()
      hc <- hclust(dist(m,method=input$dmeth), method=input$meth)
      max_h <- max(hc$height)
      sliderInput("cutval","Choose cut value:", value = max_h/2,
                  min = 0, max = max_h, step = 1)
    }
  })
  
  #create the button Select/Deselect pattern groups
  output$select_all <- renderUI({
    if(!is.null(input$datafile))
      actionButton("select_all", "Select/Deselect all")
  })
  
  
  #main page: show some instructions ....
    output$main <- renderUI({
      str1 <- paste("Main functions: ")
      str2 <- paste("1. Clustering  ")
      str3 <- paste("2. Visualization  ")
      str4 <- paste("3. Detail  ")
      HTML(paste(str1, str2, str3, str4, sep = '<br/>'))
    })
    
    #creat hierarchical cluster
    output$cluster <- renderPlot({
      m <- read_pattern() #get pattern set
      #creat hierarchical cluster according to the selected clustering method and distance measure
      hc <- hclust(dist(m,method=input$dmeth), method=input$meth) 
      plot(hc, xlab=paste(input$dmeth, "distance; ", input$meth, "clustering"))
      max_cluster_height <<- max(hc$height)
      #select the cutval to cut the tree
      abline(h=input$cutval, lty=2, col="red") 
    })
   
    #creat pattern groups
    create_mp_new <- reactive({
      #load data
      if(!is.null(input$cutval))
      {
        m <- read_pattern()
        #cut hierarchical based on cutval to get pattern groups
        hc = hclust(dist(m,method=input$dmeth), method=input$meth)
        clusterCut <- cutree(hc, h=input$cutval) #cut based on height value
        #clusterCut <- cutree(hc,k=10) #cut based on number of clusters
        
        x <- cbind(data.frame(m),clusterCut)
        da <- m
        
        ###create pattern groups by union, intersection, or majority
        # create pattern groups
        lmda <- list() #list contains pattern groups
        mt <- input$mp_method #get method from radioButtons: union, intersection, majority
        switch(mt,
               "union"={#print("Union")
                 for(i in 1:max(clusterCut)){ #for each sub-cluster
                   sp <- subset(x,clusterCut==i) #get all rows in a sub-cluster i
                   a <- c()
                   for(j in 1:nrow(sp)){ #for each row in sub-cluster
                     b <- da[(as.integer(row.names(sp[j,]))),] #get a row for da to b
                     a <- union(a,b) #union a and b  
                   }
                   lmda[[i]] <- unique(a)
                 }
               },
               "intersection" ={#print("Intersection")
                 for(i in 1:max(clusterCut)){ #for each sub-cluster
                   sp <- subset(x,clusterCut==i) #get all rows in a sub-cluster i
                   a <- c(da[(as.integer(row.names(sp[1,]))),])
                   if(nrow(sp)>1)
                   {
                     for(j in 2:nrow(sp)){ #for each row in sub-cluster
                       b <- da[(as.integer(row.names(sp[j,]))),] #set a row from da to b
                       a <- intersect(a,b) #intersection a and b
                     }  
                   }
                   lmda[[i]] <- unique(a)
                 }
               },
               "majority" = {#print("Majority")
                 for(i in 1:max(clusterCut)){ #for each sub-cluster
                   sp <- subset(x,clusterCut==i) #get all rows in a sub-cluster i
                   
                   a <- c()
                   for(j in 1:nrow(sp)){ #for each row in sub-cluster
                     b <- da[(as.integer(row.names(sp[j,]))),] #set a row from da to b
                     a <- c(a,b)
                   }
                   ## keep SNP presenting at least fre times
                   aa <- c()
                   t <- 1 
                   fre <- input$fre*nrow(sp)
                   for(k in 1:length(a))
                   {
                     if(length(grep(a[k],a))>fre && !is.na(a[k]))
                     {
                       aa[t] <- a[k]
                       t <- t+1
                     }
                   }
                   lmda[[i]] <- unique(aa)
                 }
                 
               }
        )
        
        lmda #return list of pattern groups
        
      } #end if
      
    
    })
    
    #create checkboxgroups to select pattern groups
    output$mp_checkgroup <- renderUI({
      if(!is.null(input$datafile) && !is.null(input$cutval))
      {
        m <- create_mp_new() #create pattern groups
        mp_label <- c()      #init checkbox label
        mp_id <- c()         #init checkbox ID
        #create label and ID for each pattern group
        for(i in 1:length(m))
        {
          mp_label[[i]] <- paste("Pattern group ",i)
          mp_id[[i]] <- i
        }
        #create the checkboxgroups
        checkboxGroupInput("mp_groups", "Choose pattern groups:",
                           choiceNames = mp_label,
                           choiceValues = mp_id,
                           selected = NULL)
      }
    })
  
    #update Selec/Deselect button
    observe({
      if(!is.null(input$select_all))
      {
        if (input$select_all > 0) {
          m <- create_mp_new()
          mp_label <- c()
          mp_id <- c()
          for(i in 1:length(m))
          {
            mp_label[[i]] <- paste("Pattern group ",i)
            mp_id[[i]] <- i
          }
          if(input$select_all %%2 == 0)
          {
            updateCheckboxGroupInput(session=session, inputId="mp_groups",
                                     choiceNames = mp_label,
                                     choiceValues = mp_id,
                                     selected=NULL)
          }
          else
          {
            updateCheckboxGroupInput(session=session, inputId="mp_groups",
                                     choiceNames = mp_label,
                                     choiceValues = mp_id,
                                     selected=mp_id)
          }
        }
        
      }
      
    })

    #create plot for pattern groups
    output$meta_plot <- renderPlot({
      
      #read chromosome
      ch <- read_chr()           #read chromosome detail
      max_chr_len <- max(ch[,2]) #get maximal length of chromosomes
      ch_name <- ch[,1]          #get chromosome names
      ch_len <- ch[,2]           #get chromosome lengths
      ch_centromere <- ch[,3]    #get chromosome centromere
      nb_chr <- nrow(ch)         #get number of chromosomes
      wid <- max_chr_len/(nb_chr+nb_chr/2) #calculate chromosome width
      gap <- wid/2               #the gap (space) between 2 chromosomes
      x <- c(0,max_chr_len)      #x axis and y axis
      y <- c(0,max_chr_len)
      
      #creat ytick and ylab = ch_name
      ytick <- c()
      ylabs <- c()
      for(i in 1:nb_chr){
        ytick[i] <-max_chr_len - (i-1)*(wid+gap)
        ylabs[i] <- paste(ch_name[i])  
      }
      
      #read QTL 
      qtl <- read_QTL()
      qtl_ch <- qtl[,1]
      qtl_name <- qtl[,2]
      qtl_start <- qtl[,3]
      qtl_stop <- qtl[,4]
      
      
      #draw an empty plot      
      plot(x,y,main="Overview",type="n",xaxt="n",yaxt="n",xlab="",ylab="")
      axis(2, at=ytick,labels=ylabs, col.axis="blue", las=2)
      axis(1, col.axis="blue", las=0)
      
      #draw empty chromosomes
      for(i in 1:length(ch_len))
      {
        r <- 1000000 #length of round conner
        y_mid <- max_chr_len - (i-1)*(wid+gap)
        
        left_x_mid <- ch_centromere[i]/2
        left_size <- left_x_mid-r
        
        right_x_mid <- ch_centromere[i] + (ch_len[i]-ch_centromere[i])/2
        right_size <- right_x_mid - ch_centromere[i]-r
        
        #draw left centromere
        roundrect(mid=c(left_x_mid,y_mid), left_size, wid/2, r , dr = 0.01,
                 col = "gray95", lcol = "black", lwd = 1, angle = 0)
      
        #draw right centromere
        roundrect(mid=c(right_x_mid,y_mid), right_size, wid/2, r, dr = 0.01,
                 col = "gray95", lcol = "black", lwd = 1, angle = 0)
        
        #draw QTLs for RapeSeed only
        if(!is.null(read_QTL))
        {
          for(j in 1:length(qtl_ch))
          {
            #ch_name <- factor(ch_name, levels = levels(qtl_ch))
            
            qtl_ch <- factor(qtl_ch, levels = levels(ch_name))
            
            if(qtl_ch[j] == ch_name[i])
              segments(qtl_start[j], y_mid-1.2*gap, qtl_stop[j],y_mid-1.2*gap,col="red" )
          }
        }
        #############################
      }
      
      #draw SNPs of pattern groups on chromosomes
      if(!is.null(input$mp_groups))
      {
        pts <- create_mp_new() #get pattern set = pts
        index <- read_index() #read index file
        nb_mp <- 1
        #draw legend label
        text(wid,0,"PG colors: ",col="blue")
        draw_post <- 2*wid
        for(k in input$mp_groups) #for each selected pattern group
        {
          i <- strtoi(k)
          #if(!is.null(pts[[i]]))
          {
            p_color <- i+1 #color of pattern group i^th
            #draw legend color ids
            text(draw_post+i*wid,0,k,col=p_color)
            nb_mp <- nb_mp+1
            #pt <- unique(pts[[i]]) #get pattern group i^th
            pt <- pts[[i]]
            for(j in pt){ #for each item (SNP) in pattern group i^th
              it <- strtoi(j) #get item j^th
              if(!is.na(it))
              {
                snp_detail <- index[it+1,]#get snp information detail
                snp_chr <- snp_detail[,4] #get chromosome name
                snp_pos <- snp_detail[,5] #get chromosome position
                if(!is.na(snp_chr))
                {
                  for(l in 1:nb_chr) #draw snp on correspond chromosome
                  {
                    if(as.character(snp_chr) == as.character(ch_name[l]))
                      segments(snp_pos, ytick[l]-wid/2, snp_pos, ytick[l]+wid/2, col= p_color)
                  }
                }
              }
            }#end if null
            
          }
          
        }
      }#endif
    }) #end plot
    
    #Show position on chromosome when mouse click on pattern group plot
    output$meta_pos_info <- renderText({
      if(is.null(input$meta_plot_click))
        paste0("position : ")
      else
        paste0("position : ", round(input$meta_plot_click$x,0)/1000000," Mbp")
    })
    
    #Show position of chromosome when mouse click on zoom-meta-plot
    output$zoom_meta_pos_info <- renderText({
      if(is.null(input$zoom_meta_plot_click))
        paste0("position : ")
      else
        paste0("position : ", round(input$zoom_meta_plot_click$x,0)/1000000," Mbp")
    })

    #get ranges (selected region) to zoom
    meta_ranges <- reactiveValues(x = NULL, y = NULL)
    observeEvent(input$meta_plot_bdlclick, {
      brush <- input$meta_brush
      if (!is.null(brush))
        {
        meta_ranges$x <- c(brush$xmin, brush$xmax)
        meta_ranges$y <- c(brush$ymin, brush$ymax)
        }
      else
        {
        meta_ranges$x <- NULL
        meta_ranges$y <- NULL
        }
    })

    #zoom the selected region
    output$zoom_meta_plot <- renderPlot({
      
      #get selected region
      xranges <- meta_ranges$x
      yranges <- meta_ranges$y
      
      # indicate chromosome base on yranges
      #if(!is.null(input$meta_plot_bdlclick))
      {
        ch <- read_chr()           #read chromosome detail
        max_chr_len <- max(ch[,2]) #get maximal length of chromosomes
        ch_name <- ch[,1]
        ch_len <- ch[,2]
        ch_centromere <- ch[,3]
        nb_chr <- nrow(ch)         #get number of chromosomes
        wid <- max_chr_len/(nb_chr+nb_chr/2)       #chromosome width
        gap <- wid/2               #gap between 2 chromosomes
        mp_pos <- max_chr_len/1000 #position to draw pattern group id in chromosome
        x <- c(0,max_chr_len)
        y <- c(0,max_chr_len)
        
        #creat ytick and ylab
        ytick <- c()
        ylabs <- c()
        for(i in 1:nb_chr){
          ytick[i] <-max_chr_len - (i-1)*(wid+gap)
          ylabs[i] <- paste(ch_name[i])  
        }
        
        #draw an empty plot      
        plot(x,y,main="Zoom",type="n",xaxt="n",yaxt="n",xlab="",ylab="",
             xlim=meta_ranges$x,ylim=meta_ranges$y)
        
        axis(2, at=ytick,labels=ylabs, col.axis="blue", las=2)
        axis(1, col.axis="blue", las=0)
        
        #read QTL 
        qtl <- read_QTL()
        qtl_ch <- qtl[,1]
        qtl_name <- qtl[,2]
        qtl_start <- qtl[,3]
        qtl_stop <- qtl[,4]
        
        #draw empty chromosomes
        for(i in 1:length(ch_len))
        {
          r <- 1000000 #length of round conner
          y_mid <- max_chr_len - (i-1)*(wid+gap)
          wid_des <- nb_chr #decrease width size

          left_x_mid <- ch_centromere[i]/2
          left_size <- left_x_mid-r

          right_x_mid <- ch_centromere[i] + (ch_len[i]-ch_centromere[i])/2
          right_size <- right_x_mid - ch_centromere[i]-r
          #draw left centromere
          roundrect(mid=c(left_x_mid,y_mid), left_size, wid/wid_des, r , dr = 0.01,
                    col = "gray95", lcol = "black", lwd = 1, angle = 0)
          #draw right centromere
          roundrect(mid=c(right_x_mid,y_mid), right_size, wid/wid_des, r, dr = 0.01,
                    col = "gray95", lcol = "black", lwd = 1, angle = 0)
          
          #draw QTLs for RapeSeed only
          if(!is.null(read_QTL))
          {
            for(j in 1:length(qtl_ch))
            {
              #ch_name <- factor(ch_name, levels = levels(qtl_ch))
              qtl_ch <- factor(qtl_ch, levels = levels(ch_name))
              
              if(qtl_ch[j] == ch_name[i])
              {
                segments(qtl_start[j], y_mid-1.3*gap, qtl_stop[j],y_mid-1.3*gap,col="red" )
                text(qtl_start[j]+(qtl_stop[j]-qtl_start[j])/2,y_mid-1.4*gap,qtl_name[j])

              }

            }
          }
          ###########
          
        }
 
        # #draw gene information: time-consuming because the gene file is very large
        # if(!is.null(input$meta_plot_bdlclick))
        # {
        #   #read gene info
        #   gene_info <- read_gene()
        #   gen_ch <- gene_info[,1]
        #   gen_start <-gene_info[,2]
        #   gen_stop <-gene_info[,3]
        #   gen_name <- gene_info[,4]
        # 
        #   curr_pos <- max_chr_len
        #   #get selected chromosome
        #   tt <- 1
        #   for(i in 1:nrow(ch))
        #   {
        #     if(yranges[2] < curr_pos)
        #     {
        #       curr_pos <- curr_pos - wid*3/2
        #       tt <- tt+1
        #     }
        #     else
        #     {
        #       break()
        #     }
        #   }
        #   #draw gene on selected chromosome
        #   for (g in 1:length(gen_ch) )
        #   {
        #     ch_name <- factor(ch_name, levels = levels(gen_ch))
        #     if(ch_name[tt] == gen_ch[g])
        #     {
        #       #segments(strtoi(gen_start[g]),strtoi(ch_len[tt]),strtoi(gen_stop[g]),strtoi(ch_len[tt]),col = "red")
        #       segments(gen_start[g],curr_pos-gap,gen_stop[g],curr_pos-gap,col = "red")
        #       text(strtoi(gen_start[g]),curr_pos-(3/2)*gap,gen_name[g],srt=90)
        #       #print(gen_name[g])
        #     }
        #   }
        # }
        # 
        # 
        
        #draw selected SNPs on chromosomes
        if(!is.null(input$mp_groups))
        {
          pts <- create_mp_new() #get pattern group set = pts
          index <- read_index() #read index file
          nb_mp <- 1
          
          #read gene info
          gene_info <- read_gene()
          gen_ch <- gene_info[,1]
          gen_start <-gene_info[,2]
          gen_stop <-gene_info[,3]
          gen_name <- gene_info[,4]

          for(k in input$mp_groups) #for each selected pattern group
          {
            i <- strtoi(k)
            #if(!is.null(pts[[i]]))
            {
              p_color <- i+1 #color of pattern group i^th
              #pt <- unique(pts[[i]]) #get pattern group i^th
              pt <- pts[[i]] #get pattern group i^th
              for(j in pt){ #for each item in pattern group i^th
                it <- strtoi(j) #get item j^th
                if(!is.na(it))
                {
                  snp_detail <- index[it+1,]#get snp information detail
                  snp_name <- snp_detail[,2]
                  snp_allel <- snp_detail[,3]
                  snp_chr <- snp_detail[,4] #get chromosome name
                  snp_pos <- snp_detail[,5] #get chromosome position
                  if(!is.na(snp_chr))
                  {
                    for(l in 1:nb_chr) #draw snp on correspond chromosome
                    {
                      if(as.character(snp_chr) == as.character(ch_name[l]))
                      {
                        segments(snp_pos, ytick[l]-wid/wid_des, snp_pos, ytick[l]+wid/wid_des, col= p_color)
                        text(snp_pos,ytick[l]-wid/wid_des-gap/2,paste0(snp_name,"_",snp_allel) ,srt=90,adj = 0.1)
                        text(snp_pos,ytick[l]+wid/wid_des+(nb_mp+1-1)*mp_pos,k,col=p_color)
                        
                        ###find gene containing snp_pos ??? and draw that gene
                          for (g in 1:length(gen_ch) )
                          {
                            #ch_name <- factor(ch_name, levels = levels(gen_ch))
                            #if(ch_name[tt] == gen_ch[g])
                            if( (snp_pos >= gen_start[g])&& (snp_pos <= gen_stop[g] ))
                            {
                              #segments(strtoi(gen_start[g]),strtoi(ch_len[tt]),strtoi(gen_stop[g]),strtoi(ch_len[tt]),col = "red")
                              segments(gen_start[g],ytick[l],gen_stop[g],ytick[l],col = "red")
                              text(strtoi(gen_start[g]),ytick[l]+wid/2,gen_name[g],srt=90)
                              break()
                            }
                          }
                        ###################
                        
                        
                      }
                    }
                  }
                }
              }#endif null
              
            }
            nb_mp = nb_mp+1
          }
        } #endif null
      }
    }) #end meta plot zoom
    
    #view pattern group in a detailed
    output$meta_pattern_detail <- renderDataTable({
      
      pts <- create_mp_new() #get pattern group sets
      index <- read_index() #get index data
      new_pts <- list() #new_pts for snp_name 
      #creat a table for pattern group with snp_name-pos-chr  
      if(!is.null(input$mp_groups))
      {
        n <- 1 #count number of pattern groups
        for(k in input$mp_groups)
        {
          i <- strtoi(k)
            mp <- pts[[i]] #get pattern group i^th
            if(NCOL(na.omit(mp))==0) #check if mp is empty??
            {
              new_pts[[n]] <- c("NULL")
            }
            else
            {
              new_mp <- c()
            t <- 1 #count index
            for(j in mp){ #for each item in pattern group i^th
              it <- strtoi(j) #get item j^th
              snp_detail <- index[it+1,] #get chr and pos of a SNP at index j^th
              snp_name <- snp_detail[,2]
              snp_allel <- snp_detail[,3]
              snp_chr <- snp_detail[,4]
              snp_pos <- snp_detail[,5]
              snp <- paste(snp_chr,"_",snp_pos,"_",snp_name,"_",snp_allel)
              new_mp[t] <- as.character(snp)
              t <- t+1
            }
            new_pts[[n]] <- new_mp
            }
            n <- n+1
        }
      }

      #transform list of list to data frame
      max_col <- sapply(new_pts,length)

      res <- as.data.frame(do.call(rbind,lapply(new_pts,'length<-',max(max_col))))
      #transpose data frame and print
      meta_pattern <- t(res)
      
      #create column name
      c_name <- c()
      k=1
      for(i in input$mp_groups)
      {
        c_name[k] <- paste0("Pattern group ",i)
        k=k+1
      }
      colnames(meta_pattern) <- c_name
      
      data.frame(meta_pattern)
    })
    

    

    #############END###########
}

