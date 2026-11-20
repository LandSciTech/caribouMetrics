#' Plot estimated values from Bayesian model results
#' 
#' 
#'
#' @inheritParams trajectoriesFromBayesian
#' @param metrics character. Metrics to plot, can be Recruitment or Survival or both
#' @param plotMissing logical. Should years with missing input data have interpolated results plotted?
#' @param y_labs character. Names for the y axis.
#'
#' @returns a ggplot object
#' @export
#'
#' @examples
#' surv_data <- bboudata::bbousurv_a %>% filter(Year > 2010)
#' recruit_data <- bboudata::bbourecruit_a %>% filter(Year > 2010)
#' bbouInformative <- estimateBayesianRates(surv_data, recruit_data,
#'                                          return_mcmc = TRUE, 
#'                                          # only set to speed up examples
#'                                          niters = 20) 
#' 
#' plotMeans(bbouInformative)
plotMeans <- function(bayesianResults,
                      metrics = c("Recruitment", "Survival"),
                      plotMissing = FALSE,
                      y_labs = paste0("Estimated ", metrics),
                      i18n = NULL,
                      sess = shiny::getDefaultReactiveDomain()){
  if(is.null(i18n)){
    i18n <- list(t = function(x, session)paste0(x),
                 get_languages = function(x)"en")
  } 
  
  if(inherits(bayesianResults, "data.frame")){
    stop("plotMeans only works when return_mcmc is true")
  }
  
  dat <- bayesianResults$parList
  plts <- purrr::map2(metrics, y_labs, function(metric, y_lab){
    if(dat$type == "beta"){
      dat_x <- dat[[metric]]
      
      if(is.null(dat_x)) stop("metric ", metric, " is not an element of bayesianResults$parList")
      
      # remove rows for monthly if not actually monthly
      dat_x <- dat_x %>% distinct(mean, sd, lower, upper, Annual, PopulationName, .keep_all = TRUE)
      if(!plotMissing){
        if(metric == "Survival"){
          wch_metric <- c("surv_fit", "Mortalities")
        } else {
          wch_metric <- c("recruit_fit", "Cows")
        }
        
        dat_x <- dat_x %>% 
          semi_join(bayesianResults[[wch_metric[1]]]$data %>% filter(!is.na(.data[[wch_metric[2]]])),
                    by = join_by(Annual, PopulationName))
      }
      plt <- dat_x %>%
        ggplot2::ggplot(ggplot2::aes(x = Year, y = mean, ymin = lower, ymax = upper))+
        ggplot2::geom_pointrange()+
        ggplot2::facet_wrap(~PopulationName)+
        ggplot2::labs(y = y_lab)
      
    } else {
      
      if(metric == "Survival"){
        plt <- bboutools::bb_plot_year_survival(bayesianResults$surv_fit)
        
        if(!plotMissing){
          plt$data <- plt$data %>%
            left_join(bayesianResults$surv_fit$data %>%
                        mutate(PopulationName, Annual = as.character(Annual) %>%
                                 as.numeric(), has_data = 1, .keep = "used"),
                      by = c("PopulationName", CaribouYear = "Annual")) %>%
            mutate(has_data = replace_na(has_data, 0) %>% as.character()) %>%
            filter(has_data == "1")
        }
        
      } else if (metric == "Recruitment"){
        plt <- bboutools::bb_plot_year_calf_cow_ratio(bayesianResults$recruit_fit)
        
        if(!plotMissing){
          plt$data <- plt$data %>%
            left_join(bayesianResults$recruit_fit$data %>%
                        mutate(PopulationName, Annual = as.character(Annual) %>%
                                 as.numeric(), has_data = 1, .keep = "used"),
                      by = c("PopulationName", CaribouYear = "Annual")) %>%
            mutate(has_data = replace_na(has_data, 0) %>% as.factor()) %>%
            filter(has_data == "1")
          
        }
      } else {
        stop("metric ", metric, " is not supported. ",
             "Only Survival and Recruitment are supported")
      }
      
    }
    if(metric == "Survival"){
      plt <- plt +
        ggplot2::labs(x = i18n$t(session = sess, "Year"))+
        ggplot2::scale_x_continuous(breaks = scales::extended_breaks(5, Q = 1:5,
                                                                     w = c(0.25, 0.2, 0.1, 0.5)))+
        ggplot2::scale_y_continuous(i18n$t(session = sess, y_lab),
                                    labels = scales::percent)+
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                       legend.position = "top")
    }else if (metric == "Recruitment"){
      plt <- plt +
        ggplot2::scale_y_continuous(transform = scales::new_transform(
          "hundred", transform = function(x){x*100}, inverse = function(x){x/100},
          format = scales::label_number(scale = 100)
        ))+
        ggplot2::scale_x_continuous(breaks = scales::extended_breaks(5, Q = 1:5,
                                                                     w = c(0.25, 0.2, 0.1, 0.5)))+
        ggplot2::labs(x = i18n$t(session = sess, "Year"), 
                      y = i18n$t(session = sess, y_lab))+
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                       legend.position = "top")
    }
    
    plt
    
  })
  return(plts)  
  
}

