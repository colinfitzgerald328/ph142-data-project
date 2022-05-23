## Part 1

1.  \[2 marks\] What is the problem your are addressing with these data?
    State the question you are trying to answer and let us know what
    type of question this is in terms of the PPDAC framework.

The problem we are addressing with this data is the estimation of
obesity levels based on eating habits and physical condition.

We are trying to answer the question of whether there are specific
factors within eating habits and physical conditions that leads to a
higher levels of obesity. Some of these factors include frequent
consumption of high caloric food, number of main meals, consumption of
water daily, physical activity frequency, and transportation used. This
question is a causative/etiologic question:

Are there potential associations or trends within the eating habits and
physical conditions of individuals of Mexico, Peru and Columbia that may
have lead to higher levels of obesity in the respective countries?

1.  \[2 marks\] What is the target population for your project? Why was
    this target chosen i.e., what was your rationale for wanting to
    answer this question in this specific population?

The target population for this project are individuals from Mexico,
Peru, and Columbia between the ages of 14 and 61 with diverse eating
habits and physical conditions.

The rationale for wanting to answer this question in this specific
population is that we want to study potential community or cultural
factors that may have led to high obesity levels and be able to make a
generalization of risk factors for these specific countries.

There are people with all different factors and exposures that can help
us identify key factors that lead to an increase in obesity levels.
Specifically, these countries were chosen as they have high rates of
obesity and are currently showing no steps at lowering them. This may
provide insight on the factors that increase obesity levels and may show
us areas that can be focused on to lower these levels.

1.  \[2 marks\] What is the sampling frame used to collect the data you
    are using? Describe why you think this sampling strategy is
    appropriate for your question. To what group(s) would you feel
    comfortable generalizing the findings of your study and why.

This data comes from an estimation of obesity levels in people from the
countries of Mexico, Peru and Colombia, with ages between 14 and 61 and
diverse eating habits and physical condition.

Data was collected using a web platform with a survey where anonymous
users answered each question, then the information was processed
obtaining 17 attributes and 2111 records. The attributes that were
collected were related to eating habits, physical conditions, gender,
age, height and weight. This sampling strategy is appropriate for our
question given we are interested in looking at and getting more
information on potential community or culture driven factors that lead
to a higher levels of obesity in each country.

We would feel comfortable generalizing the findings of our study to
Mexico, Peru, and Colombia’s neighboring countries in Central and South
America given that we may find potential trends and associations in the
eating and physical habits across the three countries.

1.  \[2 marks\] Write a brief description (1-4 sentences) of the source
    and contents of your dataset. Provide a URL to the original data
    source if applicable. If not (e.g., the data came from your
    internship), provide 1-2 sentences saying where the data came from.
    If you completed a web form to access the data and selected a
    subset, describe these steps (including any options you selected)
    and the date you accessed the data.

<https://archive.ics.uci.edu/ml/datasets/Estimation+of+obesity+levels+based+on+eating+habits+and+physical+condition+>.

The data came from the UC Irvine Machine Learning Repository and
contains information such as Frequent consumption of high caloric food
(FAVC), Frequency of consumption of vegetables (FCVC), Physical activity
frequency (FAF), Time using technology devices (TUE), Transportation
used (MTRANS), and other variables as such. We didn’t have to complete a
web form to access the data and did not need to subset it either. There
were links to download both the Data Folder and the Data Set Description
at the top of the website. We accessed the data on Monday February 21,
2022.

1.  \[1 mark\] Write code below to import your data into R. Assign your
    dataset to an object.

<!-- -->

    obesity_data <- read_csv("ObesityDataSet_raw_and_data_sinthetic.csv")

1.  \[3 marks\] Use code in R to answer the following questions:

<!-- -->

1.  What are the dimensions of the dataset?

<!-- -->

    dim(obesity_data)

    ## [1] 2111   17

The dimensions of the dataset are 2111 rows by 17 columns.

1.  Provide a list of variable names.

<!-- -->

    names(obesity_data)

    ##  [1] "Gender"                         "Age"                           
    ##  [3] "Height"                         "Weight"                        
    ##  [5] "family_history_with_overweight" "FAVC"                          
    ##  [7] "FCVC"                           "NCP"                           
    ##  [9] "CAEC"                           "SMOKE"                         
    ## [11] "CH2O"                           "SCC"                           
    ## [13] "FAF"                            "TUE"                           
    ## [15] "CALC"                           "MTRANS"                        
    ## [17] "NObeyesdad"

1.  Print the first six rows of the dataset.

<!-- -->

    head(obesity_data)

    ## # A tibble: 6 × 17
    ##   Gender   Age Height Weight family_history_with_… FAVC   FCVC   NCP CAEC  SMOKE
    ##   <chr>  <dbl>  <dbl>  <dbl> <chr>                 <chr> <dbl> <dbl> <chr> <chr>
    ## 1 Female    21   1.62   64   yes                   no        2     3 Some… no   
    ## 2 Female    21   1.52   56   yes                   no        3     3 Some… yes  
    ## 3 Male      23   1.8    77   yes                   no        2     3 Some… no   
    ## 4 Male      27   1.8    87   no                    no        3     3 Some… no   
    ## 5 Male      22   1.78   89.8 no                    no        2     1 Some… no   
    ## 6 Male      29   1.62   53   no                    yes       2     3 Some… no   
    ## # … with 7 more variables: CH2O <dbl>, SCC <chr>, FAF <dbl>, TUE <dbl>,
    ## #   CALC <chr>, MTRANS <chr>, NObeyesdad <chr>

1.  \[4 marks\] Use the data to demonstrate a statistical concept from
    Part I of the course. Describe the concept that you are
    demonstrating and interpret the findings. This should be a
    combination of code and written explanation.

The overarching statistical concept that our group wanted to demostrated
in the linear association between two quantitative variables from our
dataset. The two variables we did this with were height of participants
(in meters) and weight of participants (in kilograms). The method we
used to do so was performing linear regression of the dependent variable
in this case of weight on the explanatory variable of height. We first
made a scatter plot through ggplot2 with a best fit regression line
include from the first code chunk below to visualize the data and get a
sense of what the regression coefficients should look like. We then
performed linear regression shown in the second code chunk below and
received a readout for the regression coefficients as -134.6 for the
intercept coefficient and 130.0 for the coefficient of height on weight.
In context, the intercept coefficient indicates that a person with a
height of 0.0 meters should weigh -134.6 kg. Obviously this doesn’t make
sense in context because people cannot have negative weights and this
readout comes from the narrowness of our sampling data. Only people
between roughly 1.4 meters and 2.0 meters in height are included in this
dataframe, and because this data is looking at primarily obese
individuals, the data most likely has a somewhat unrepresentative height
to weight ratio when performing regression as compared to the general
population. In context, the readout for the height coefficient suggests
that a 1 meter increase in the height of a participant will mean that
there will be a 130 kg increase in the weight of the individual. This
coefficient makes more sense in context as there tends to be a direct
relationship between height and weight in people, but the magnitude of
this coefficient is likely biased to the selection of the sampling frame
as well.

    obesity_hw_scatter = ggplot(data = obesity_data, aes(x = Height, y = Weight)) + 
      geom_point() +
      geom_smooth(method = lm) + 
      labs(title = 'Obesity Data Heights and Weights Plotted')
    obesity_hw_scatter

    ## `geom_smooth()` using formula 'y ~ x'

![](Data-Project-Part-2_files/figure-markdown_strict/unnamed-chunk-5-1.png)

    obesity_hw_regression <- lm(Weight ~ Height, data = obesity_data)
    obesity_hw_regression

    ## 
    ## Call:
    ## lm(formula = Weight ~ Height, data = obesity_data)
    ## 
    ## Coefficients:
    ## (Intercept)       Height  
    ##      -134.6        130.0

## Part 2

1.  \[2 marks\] Describe a quantity you will estimate as an outcome in
    your problem using probability notation. Are you planning to
    calculate marginal probabilities? Conditional probabilities?

<!-- -->

    library(janitor)

    ## Warning in system("timedatectl", intern = TRUE): running command 'timedatectl'
    ## had status 1

    ## 
    ## Attaching package: 'janitor'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     chisq.test, fisher.test

    obesity_data$CAEC <- factor(obesity_data$CAEC, levels=c('no', 'Sometimes', 'Frequently', 'Always'))
    tabyl(obesity_data, Gender, CAEC) %>% adorn_totals(where = c("row", "col")) %>% 
      adorn_percentages("row") %>% adorn_pct_formatting() %>% adorn_ns() %>% knitr::kable()

<table>
<colgroup>
<col style="width: 10%" />
<col style="width: 15%" />
<col style="width: 19%" />
<col style="width: 18%" />
<col style="width: 15%" />
<col style="width: 21%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Gender</th>
<th style="text-align: left;">no</th>
<th style="text-align: left;">Sometimes</th>
<th style="text-align: left;">Frequently</th>
<th style="text-align: left;">Always</th>
<th style="text-align: left;">Total</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Female</td>
<td style="text-align: left;">1.4% (15)</td>
<td style="text-align: left;">80.9% (844)</td>
<td style="text-align: left;">15.4% (161)</td>
<td style="text-align: left;">2.2% (23)</td>
<td style="text-align: left;">100.0% (1043)</td>
</tr>
<tr class="even">
<td style="text-align: left;">Male</td>
<td style="text-align: left;">3.4% (36)</td>
<td style="text-align: left;">86.2% (921)</td>
<td style="text-align: left;">7.6% (81)</td>
<td style="text-align: left;">2.8% (30)</td>
<td style="text-align: left;">100.0% (1068)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Total</td>
<td style="text-align: left;">2.4% (51)</td>
<td style="text-align: left;">83.6% (1765)</td>
<td style="text-align: left;">11.5% (242)</td>
<td style="text-align: left;">2.5% (53)</td>
<td style="text-align: left;">100.0% (2111)</td>
</tr>
</tbody>
</table>

We would like to calculate the conditional probabilities of eating
habits based on gender. We are curious to see how eating habits between
meals vary by gender. The above table creates those conditional
probabilities for us. The quantities we will estimate as an outcome, are
categorical variables, (no, sometimes, frequently, always), based on
gender. Below we report those calculated probabilities. These are
estimations for this data set only.

*P*(does not eat between meals | Female)  = 15/1043 = 1.4%

*P*(sometimes eats between meals | Female)  = 844/1043 = 80.9%

*P*(frequently eats between meals | Female)  = 161/1043 = 15.4%

*P*(always eats between meals | Female)  = 23/1043 = 2.2%

*P*(does not eat between meals | Male)  = 36/1068 = 3.4%

*P*(sometimes eats between meals | Male)  = 921/1068 = 86.2%

*P*(frequently eats between meals | Male)  = 81/1068 = 7.6%

*P*(always eats between meals | Male)  = 30/1068 = 2.8%

The above probabilities are the probability of each eating habit
conditioned on gender.

1.  \[3 marks\] Describe the type of theoretical distribution that is
    relevant for your data.

<!-- -->

1.  What type of variable(s) are you investigating (continuous,
    categorical, ordinal, etc)?

There are 17 attributes in our data set and all of the variables were
categorical except three: Weight, Height, and Age. Those three variables
were numerical and continuous since we can assume an infinite number of
real values within a given interval. It is worth noting that our data
set also includes 5 binary attributes since each observation either
falls into YES or NO categories or in the case of Gender ( MALE or
FEMALE), there are a fixed number of observations, each observation is
independent, and the probability of each observation is 50%:

Gender Frequent consumption of high caloric food (FAVC) Calories
consumption monitoring (SCC) Smoker Familial Overweight Connection

Some other variables of interest are Time using technology devices
(TUE), Number of main meals (NCP), Consumption of food between meals
(CAEC), which are all categorical and ordinal variables in our data set.
This is the case since their values are defined by an order relation
between the different categories ie. sometimes, frequently, and always.

1.  What theoretical distribution that we have talked about would
    potentially be appropriate to use with these data (Normal, Binomial,
    Poisson. . . )

If we were to investigate one of the continuous variables such as height
or age, then a normal distribution would be appropriate. However, if we
were to investigate one of the binary variables then a binomial
distribution would be more appropriate, which could be modeled from a
Bernoulli variable derived from setting ‘yes’ and ‘no’ to 0 and 1 in the
data set.

The main distributions for categorical data analysis are the binomial
and multinomial distributions. Probabilities associated with different
combinations of events can be determined based on the distributions.

1.  Why is this an appropriate model for the data you are studying?

Let’s focus on investigating height for the purpose of this question. In
order to check whether or not the normal distribution is appropriate, we
plotted a QQ plot.

    ggplot(obesity_data, aes(sample = Height)) + stat_qq() + stat_qq_line() +
    theme_minimal(base_size = 15) + 
      labs(title = "QQ Plot of Height from Obesity Data Set") + 
      xlab("Theoretical") + ylab("Sample")

![](Data-Project-Part-2_files/figure-markdown_strict/unnamed-chunk-8-1.png)

When we plot the theoretical quantiles on the Q-Q plot, we can see that
the points are most closely attached from about a height of 1.6 to 1.8
meters so a normal distribution is appropriate for those values;
however, as we go higher and lower we start to see a little bit of a
skew so a normal distribution would no longer be appropriate.

Furthermore, after stratifying height by gender, we plotted their
histograms and we can see that they follow a relatively normal
distribution for both genders.

    ggplot(obesity_data,aes(x = Height, colour=Gender)) + 
      geom_histogram(aes(y = ..count.., fill=Gender), alpha=0.5,position="identity") + 
      ggtitle("Histogram of Height by Gender from Data Set")

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Data-Project-Part-2_files/figure-markdown_strict/unnamed-chunk-9-1.png)

1.  \[4 marks\] Use the data you have to demonstrate a statistical
    concept from Part II of the course. Describe the concept that you
    are demonstrating and interpret the findings. This may include code
    in R, a visual of some kind and text interpretation

<!-- -->

    mean_height <- mean(obesity_data$Height)
    std_height <- sd(obesity_data$Height)

    men <-obesity_data %>% filter(Gender == 'Male')
    male_mean_height <- mean(men$Height)
    std_male_height <- sd(men$Height)

    women <- obesity_data %>% filter(Gender == 'Female')
    female_mean_height <- mean(women$Height)
    std_female_height <- sd(women$Height)

    CI_lower_all_data <- mean_height - (1.96 * (std_height/sqrt(nrow(obesity_data))))
    CI_upper_all_data <- mean_height + (1.96 * (std_height/sqrt(nrow(obesity_data))))

    CI_lower_male_data <- male_mean_height - (1.96 * (std_male_height/sqrt(nrow(men))))
    CI_upper_male_data <- male_mean_height + (1.96 * (std_male_height/sqrt(nrow(men))))

    CI_lower_female_data <- female_mean_height - (1.96 * (std_female_height/sqrt(nrow(women))))
    CI_upper_female_data <- female_mean_height + (1.96 * (std_female_height/sqrt(nrow(women))))

    CI_all_data <- c(CI_lower_all_data, CI_upper_all_data)
    CI_all_data

    ## [1] 1.697697 1.705658

    CI_male_data <- c(CI_lower_male_data, CI_upper_male_data)
    CI_male_data

    ## [1] 1.754362 1.763019

    CI_female_data <- c(CI_lower_female_data, CI_upper_female_data)
    CI_female_data

    ## [1] 1.638776 1.647820

    library(ggplot2)
    library(MASS)       # for fitdistr(...)

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    get.params <- function(z) with(fitdistr(z,"normal"),estimate[1:2])
    df <- aggregate(Height~Gender, obesity_data, get.params)
    df <- data.frame(Gender=df[,1],df[,2])
    x  <- with(obesity_data, seq(min(obesity_data$Height),max(obesity_data$Height),len=100))
    gg <- data.frame(Height=rep(x,nrow(df)),df)
    gg$y <- with(gg,dnorm(x,mean,sd))
    gg$y <- gg$y * aggregate(Height~Gender, obesity_data,length)$Height * diff(range(obesity_data$Height))/30

    ggplot(obesity_data,aes(x = Height, colour=Gender)) + 
      geom_histogram(aes(y = ..count.., fill=Gender), alpha=0.5,position="identity") +
      geom_line(data=gg, aes(y=y)) + 
      geom_vline(xintercept = CI_lower_female_data) + 
      geom_vline(xintercept = CI_upper_female_data) + 
      geom_vline(xintercept = CI_lower_male_data) + geom_vline(xintercept = CI_upper_male_data)

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Data-Project-Part-2_files/figure-markdown_strict/unnamed-chunk-15-1.png)

The above is a visualization of our confidence intervals, for the
stratification by gender.

The concept that we are demonstrating is finding a confidence interval.
A confidence interval displays the probability that a parameter will
fall between a pair of values around the mean. Specifically within our
data, we found the confidence interval for our entire dataset and
stratified our data by gender (male and female). Our confidence
intervals specifically work to capture the mean height of obese males,
females, and both in Mexico, Peru and Columbia. Based on our results, we
found that the interval 1.697697 - 1.705658 (meters) has a 95% success
rate in capturing within that interval the mean height µ of all obese
males and females in Mexico, Peru and Columbia. After we stratified our
data by gender we found that the interval 1.754362 - 1.763019 (meters)
has a 95% success rate in capturing within that interval the mean height
µ of obese males in Mexico, Peru and Columbia. We also found that the
interval 1.638776 - 1.647820 (meters) has a 95% success rate in
capturing within that interval the mean height µ of obese females in
Mexico, Peru and Columbia. By conducting these confidence intervals, we
found there to be a statistically significant difference in the true
mean height of all 3 intervals due to there not being any overlap in
confidence intervals.

## Part 3

1.  \[2 points\] Identify a statistical test to apply to your data (must
    be a concept we covered in part III of the course). In plain
    language, write the question you are trying to answer.

The statistical test we plan to apply to our data is a chi-square test
for independence. The question we are trying to answer is if there is an
association between obesity level and gender within individuals of our
study. This is to analyze if gender is a risk factor for obesity. We are
using 2 categories of gender (male and female) and 7 obesity levels
(insufficient weight, normal weight, obesity type 1, obesity type 2,
obesity type 3, overweight level 1, overweight level 2) to see if there
is an association between level and gender.

1.  \[2 points\] What assumptions are required by the method you chose
    in question 2? Show how you assessed whether these assumptions are
    met by your dataset.

<!-- -->

1.  An assumption that is required by the chi-squared test is to have
    performed a simple random sample to obtain the data. This data was
    collected by an online survey that may have not necessarily been
    random as we are prone to bias of response and nonresponse.
    Therefore this may not be a completely random sample; however, we
    are going to assume randomness for the purposes of this data project
    (further addressing the bias in later questions). 2) Another
    assumption is that the variables are categorical thus we made sure
    to choose 2 categorical variables - Obesity levels and Gender. With
    the obesity data, we made sure not to use weight as the determining
    factor since that is a continuous variable. Instead, we utilized the
    variable NObeyesdad, which attained obesity levels as categorical
    groups. 3) The last assumption for the chi-squared test, which our
    data met, is that the expected counts are at least 5. That can be
    clearly demonstrated in Table 3 of the Expected counts calculated
    that is presented in question 4.

1.  \[2 points\] Explain why this test is appropriate for the data you
    have and the question you are trying to answer. Use at least one
    visualization technique and include both the output and the R code
    that generated it.

A chi-square test for independence is appropriate for our data because
we are using 2 categorical variables that are able to be counted easily.
We are looking to see if there is a relationship between the two and
this test would allow us to determine if there was one. We see with our
data in the 2 way table that it is in good format for counts with a
chi-squared test and visualization with a dodge bar chart.

    ggplot(obesity_data,aes(x = NObeyesdad, colour = Gender)) + 
      geom_bar(aes(y = ..count.., fill=Gender), alpha=0.5,position = "dodge") + 
      ggtitle("Barplot of Obesity Levels by Gender from Data Set") + 
      theme(axis.text.x = element_text(angle = 45, hjust=1))

![](Data-Project-Part-2_files/figure-markdown_strict/unnamed-chunk-16-1.png)

Looking at the visualization, we can see that there appears to be a
difference in the distribution of counts of obesity level by gender,
which means a chi-squared test for independence would be appropriate to
test that presumption.

    library(janitor)
    observed_table <- matrix(c(173, 141, 156,2,323,145,103, 99,146,195,295,1,145,187), 
                             nrow = 2, ncol = 7, byrow = T)
    rownames(observed_table) <- c('Female', 'Male')
    colnames(observed_table) <- c('Insufficient Weight', 
                                  'Normal Weight', 
                                  'Obesity Type I', 
                                  'Obesity Type II',
                                  'Obesity Type III',
                                  'Overweight Level I',
                                  'Overweight Level II')
    ObservedVal_Table <- knitr::kable(observed_table, caption = "Observed Data")
    ObservedVal_Table

<table>
<caption>Observed Data</caption>
<colgroup>
<col style="width: 5%" />
<col style="width: 15%" />
<col style="width: 10%" />
<col style="width: 11%" />
<col style="width: 12%" />
<col style="width: 13%" />
<col style="width: 14%" />
<col style="width: 15%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: right;">Insufficient Weight</th>
<th style="text-align: right;">Normal Weight</th>
<th style="text-align: right;">Obesity Type I</th>
<th style="text-align: right;">Obesity Type II</th>
<th style="text-align: right;">Obesity Type III</th>
<th style="text-align: right;">Overweight Level I</th>
<th style="text-align: right;">Overweight Level II</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Female</td>
<td style="text-align: right;">173</td>
<td style="text-align: right;">141</td>
<td style="text-align: right;">156</td>
<td style="text-align: right;">2</td>
<td style="text-align: right;">323</td>
<td style="text-align: right;">145</td>
<td style="text-align: right;">103</td>
</tr>
<tr class="even">
<td style="text-align: left;">Male</td>
<td style="text-align: right;">99</td>
<td style="text-align: right;">146</td>
<td style="text-align: right;">195</td>
<td style="text-align: right;">295</td>
<td style="text-align: right;">1</td>
<td style="text-align: right;">145</td>
<td style="text-align: right;">187</td>
</tr>
</tbody>
</table>

Observed Data

    ExpectedVal_Table <- knitr::kable(chisq.test(observed_table)$expected, 
                                      caption = "Expected Data")
    ExpectedVal_Table

<table>
<caption>Expected Data</caption>
<colgroup>
<col style="width: 5%" />
<col style="width: 15%" />
<col style="width: 10%" />
<col style="width: 11%" />
<col style="width: 12%" />
<col style="width: 13%" />
<col style="width: 14%" />
<col style="width: 15%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: right;">Insufficient Weight</th>
<th style="text-align: right;">Normal Weight</th>
<th style="text-align: right;">Obesity Type I</th>
<th style="text-align: right;">Obesity Type II</th>
<th style="text-align: right;">Obesity Type III</th>
<th style="text-align: right;">Overweight Level I</th>
<th style="text-align: right;">Overweight Level II</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Female</td>
<td style="text-align: right;">134.3894</td>
<td style="text-align: right;">141.8006</td>
<td style="text-align: right;">173.4216</td>
<td style="text-align: right;">146.7414</td>
<td style="text-align: right;">160.0815</td>
<td style="text-align: right;">143.2828</td>
<td style="text-align: right;">143.2828</td>
</tr>
<tr class="even">
<td style="text-align: left;">Male</td>
<td style="text-align: right;">137.6106</td>
<td style="text-align: right;">145.1994</td>
<td style="text-align: right;">177.5784</td>
<td style="text-align: right;">150.2586</td>
<td style="text-align: right;">163.9185</td>
<td style="text-align: right;">146.7172</td>
<td style="text-align: right;">146.7172</td>
</tr>
</tbody>
</table>

Expected Data

1.  \[2 points\] Clearly state the null and alternative hypotheses for
    this test.

*H*<sub>*o*</sub>: Gender and obesity level are independent.
*H*<sub>*a*</sub>: Gender and obesity level are not independent.

1.  \[2 points\] Include the R code you used to generate your results.
    Annotate your code to help us follow your reasoning.

<!-- -->

    #perform chi-squared test in order to generate test statistic
    chisq.test(observed_table)

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  observed_table
    ## X-squared = 657.75, df = 6, p-value < 2.2e-16

    #use the chi-squared test statistic in order to 
    #calculate the p-value with the appropriate 
    #degrees of freedom (rows - 1)(columns-1) = 6
    pchisq(q=657.75, df=6, lower.tail = FALSE)

    ## [1] 8.073746e-139

1.  \[4 points\] Present your results in a clear summary. This should
    include both a text summary and a table or figure with appropriate
    labeling.

Assuming no association between gender and obesity level (NObeyesdad),
there is approximately a 0% chance of the chi-square value we calculated
or a larger one. This probability is small enough that there is evidence
in favor of the alternative hypothesis that there is a relationship
between gender and obesity level (NObeyesdad).

Below we created a visualization of the distribution of each weight
class by each gender. This table was created using the observed values
in the survey.

    library(DescTools)
    #gender_obesity_table <- tabyl(obesity_data, Gender, NObeyesdad) %>% knitr::kable()
    #ExpectedVal_Table
    #mosaicplot(gender_obesity_t, dir = c("h", "v"))
    #mosaicplot(obesity_data)
    PlotMosaic(observed_table, main = "NObeyesdad ~ Gender", 
               off = 0.05, horiz = FALSE, xlab = "NObeyesdad", y = "Gender")

![](Data-Project-Part-2_files/figure-markdown_strict/unnamed-chunk-19-1.png)

1.  \[4 points\] Interpret your findings. Include a statement about the
    evidence, your conclusions, and the generalizability of your
    findings. Our analysis and conclusions depend on the quality of our
    study design and the methods of data collection. Any missteps or
    oversights during the data collection process could potentially
    change the outcome of what we are trying to find. Consider the
    methods used to collect the data you analyzed. Was there any
    potential issue in how the participants were selected/recruited,
    retained, or assessed that may have impacted the outcome of your
    analysis/visualization? Were there any potential biases that you
    might be concerned about? Were there factors that were not measured
    or considered that you think could be important to the
    interpretation of these data?

According to the tables we generated for observed and expected counts of
gender-obesity status pairs, there were large discrepancies between the
observed and expected counts for many data pairs which is most likely
the reason we observed such an extreme chi-square test statistic and
p-value. Based on the chi-square test results, the data suggests that
there is dependence between the variables of gender and obesity status.
However, we are not particularly optimistic of the generalizability of
these findings due to the methods of data collection in order to
generate our data. The data was generated from an anonymous online
survey which could lead to two major data collection issues. The first
is related to the bias of those who chose to participate. If the survey
was targeted at viewing associations of obesity indicators to people and
their obesity status, it may have been more applicable for individuals
who fall under the obesity spectrum to participate. Additionally, even
though the survey was anonymous, people of certain groups may be more or
less likely to respond honestly to the survey given any views on their
health status.

1.  \[1 point\] Create a statement of contribution. For example, the
    American Journal of Epidemiology provides the following instructions
    to authors: “Authorship credit should be based on criteria developed
    by the International Committee for Medical Journal Editors
    (ICMJE): 1) substantial contributions to conception and design, or
    acquisition of data, or analysis and interpretation of data; 2)
    drafting the article or reviewing it and, if appropriate, revising
    it critically for important intellectual content; 3) final approval
    of the version to be published. Authors should meet all conditions.
    In addition, each author must certify that they have participated
    sufficiently in the work to believe in its overall validity and to
    take public responsibility for appropriate portions of its content.
    Author names should be listed in ScholarOne and author contributions
    should be detailed in the cover letter (e.g., “Author A designed the
    study and directed its implementation, including quality assurance
    and control. Author B helped supervise the field activities and
    designed the study’s analytic strategy. Author C helped conduct the
    literature review and prepare the Methods and the Discussion
    sections of the text.”).”

Every group member contributed relatively equally to the overarching
project given we had joint responsibilities for a lot of the questions.
Member QD had joint responsibility for question 7 in Part I of the
project, question 4 in Part II, and questions 6-8 in Part III. Member CA
had joint responsibility for questions 3 and 4 in Part I of the project,
question 3 in Part II of the project, and questions 2, 3, 4 and 5 in
Part III of the project. Member AT had joint responsibility for
questions 1 and 2 in Part I of the project, question 4 in Part II of the
project, and questions 2, 3, 4 and 5 in Part III of the project.  
Member NM had joint responsibility for questions 3 and 4 in Part I of
the project, question 3 in Part II, and questions 6-8 in Part III.
Member CF had joint responsibility for questions 5 and 6 in Part I of
the project, question 2 and 3 in part II of the project, and questions
6-8 in Part III.
