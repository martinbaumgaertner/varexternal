#' pretty_irf
#'
#' Estimates the asymptotic covariance matrix
#'
#' @param data list of var$irf objects from SVARIV
#' @param shock_names instrument names
#' @param pretty_names vector of pretty names for y-axis
#' @param cum logical. Cumulative? BEWARE: Possible mistake here
#' @param confidence_type vector of confidence types
#' @param manual_color specify colors per Instrument. Input is a named vecor with shock_names and the corresponding color (hex)
#' @param legend logical. inculde or exclude legends
#' @param title vector of headings for each column
#' @param same_scale logical value. Should all columns use the same scale?
#'
#' @return ggplot wrapper
#'
#' @examples
#'p           = 24     #Number of lags in the VAR model
#'confidence  = confidence=c(0.6,0.9,0.95)    #Confidence Level for the standard and weak-IV robust confidence set
#'NWlags      = 0;  # Newey-West lags(set it to 0 to compute heteroskedasticity robust std errors)
#'norm        = 1; # Variable used for normalization
#'scale       = 1; # Scale of the shock
#'horizons    = 20; #Number of horizons for the Impulse Response Functions(IRFs)
#(does not include the impact or horizon 0)

#load data
#'colnames(oil)<-c("a","b","c","d","year","month")
#'ydata<-oil[,1:3]
#'z<-oil[,4]
#compute IRFS
#'
#'VAR<-SVARIV(ydata,z,p,confidence,NWlags,norm,scale,horizons,instrument_name="test")
#'sh.col<-      c("#E41A1C")
#'names(sh.col)<-c("test")
#'pretty_irf(data=list(VAR$irfs),shock_names="test",pretty_names=c("a","b","c"),manual_color=sh.col,title="subheading")
#' @export
pretty_irf<-function(data,shock_names,pretty_names=NULL,cum=F,confidence_type="msw",
                     manual_color=NULL,legend=F,title=NULL,same_scale=T,
                     shock_sign="positive"
                     ){
  if(any(sapply(data, is.list))){
    number_of_shocks<-length(data)
    variable_names<-unique(as.character(data[[1]]$variable))
    lags<-unique(as.numeric(data[[1]]$horizon))-1
    conf_level<-unique(data[[1]]$confi_level)
  }else{
    number_of_shocks<-1
    variable_names<-unique(as.character(data$variable))
    lags<-unique(as.numeric(data$horizon))-1
    conf_level<-unique(data$confi_level)
  }

  if(!is.null(pretty_names)){
    variable_names_pretty<-pretty_names
  }

  if(shock_sign=="negative"){
    data<-data %>%
      dplyr::mutate(value=value*(-1))
  }

  variable_n<-length(variable_names)

    plots <- array(list(), dim = c(variable_n, number_of_shocks))
    for (i in 1:number_of_shocks){#shock
      for (j in 1:variable_n){#response
          da<-data[[i]] %>%
            dplyr::mutate(value=if(shock_sign=="negative") {value=value*(-1)
            }else{
              value=value}) %>%
            dplyr::filter(if (cum ==T) {stringr::str_detect(confi_type,"cum")
            } else {
              !stringr::str_detect(confi_type,"cum")
            })%>%
            dplyr::filter(variable==variable_names[j],
                   confi_type%in%confidence_type,
                   instrument==shock_names[i])%>%
            tidyr::pivot_wider(names_from = type, values_from = value)

          plot_temp<-ggplot2::ggplot()+
            ggplot2::geom_line(data=da,ggplot2::aes(horizon,point,group=instrument,color=instrument),size=1.2)+
            ggplot2::geom_hline(yintercept=0)+
            #ylab(variable_names_pretty[j])+
            ggplot2::xlab("")+
            ggplot2::theme_bw()+
            ggplot2::theme(plot.margin = unit(c(0,5,0,5), "mm"))+
            ylab(element_blank())

          for(a in 1:length(conf_level)){
            plot_temp<-plot_temp+
              ggplot2::geom_ribbon(data=da %>% dplyr::filter(confi_level==conf_level[a]),
                                   ggplot2::aes(horizon,ymin=lower,ymax=upper,group=instrument,fill=instrument),alpha=(0.5-a*0.1))
          }
          if(!is.null(manual_color)){#color depending on instrument
            plot_temp<-plot_temp+
              ggplot2::scale_fill_manual(values = manual_color[which(names(manual_color)==shock_names[i])])+
              ggplot2::scale_color_manual(values = manual_color[which(names(manual_color)==shock_names[i])])
          }
          if(j!=variable_n){# if not last row than no x-axis text (saved space between plots)
            plot_temp<-plot_temp+
              theme(axis.text.x=element_blank())
          }
          if(i==1){ #if first instrument include variable names
            plot_temp<-plot_temp+
              ylab(variable_names_pretty[j])
          }
          if(j==1){ # if first row include title
            plot_temp<-plot_temp+
              ggtitle(title[i])+
              theme(plot.title=element_text(size=10, vjust=-1,hjust = 0.5))
          }
          if(legend==F){
            plot_temp<-plot_temp+ theme(legend.position="none")
          }

          plots[[j,i]]<-plot_temp
      }
    }
      if(same_scale==T){# readjust y axis and remove second y axis for space
        for(i in 1:variable_n){
          maxi<-max(unlist(lapply(plots[i,], get_plot_limits)))
          mini<-min(unlist(lapply(plots[i,], get_plot_limits)))
          for(j in 1:number_of_shocks){
            plots[[i,j]]<-plots[[i,j]]+
              coord_cartesian( ylim=c(mini, maxi))
            if(j!=1){
              plots[[i,j]]<-plots[[i,j]] +
                theme(axis.title.y=element_blank(),
                      axis.text.y=element_blank())
            }
          }
        }
      }
    patchwork::wrap_plots(plots,ncol=number_of_shocks,guides="collect",byrow=F)
}
