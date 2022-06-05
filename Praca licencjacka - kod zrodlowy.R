rm(list=ls())

library(mvtnorm)
library(vars)
library(aod)

#Generowanie danych

#Zalozenia:
#Model VAR: yt = c+A*yt-1 + et

var_size = 3 #Wymiar modelu VAR
series_length = 40 #ilosc obserwacji w szeregu
series_length = series_length+100 #Dodanie 100 obserwacji dla ustabilizowania szeregu
c_const = c(0.01, 0.01, 0.01) #Stala w modelu
number_of_repetitions = 100 #Ilosc powtorzen testu

#Przedziały ufności dla asymptotycznego testu
#1%
bottom0.01 = 0.01-2*sqrt((0.01*(1-0.01))/number_of_repetitions)
top0.01 =  0.01+2*sqrt((0.01*(1-0.01))/number_of_repetitions)
#5%
bottom0.05 = 0.05-2*sqrt((0.05*(1-0.05))/number_of_repetitions)
top0.05 =  0.05+2*sqrt((0.05*(1-0.05))/number_of_repetitions)
#10%
bottom0.1 = 0.1-2*sqrt((0.1*(1-0.1))/number_of_repetitions)
top0.1 =  0.1+2*sqrt((0.1*(1-0.1))/number_of_repetitions)

#Postacie maciezry parametrow A
A1 = matrix( c(1,0,0,0,1,0,0,0,1), var_size) #Brak kointegracji
A2 = matrix( c(1,0,0.5,0,1,0.5,-0.125,0,0.5), var_size) #Dwa skointegrowane rownania
A3 = matrix( c(0.25,0,-0.75,0,1,0,-0.125,0,0.875), var_size) #Jedno skointegrowane rownanie

# ERRORS

par(mfrow=c(5,2))

#Error E1
E1 = mvrnorm(n=series_length,mu = c(0,0,0), Sigma = diag(3))
plot(E1[,1], type='l', col='red', main= "Błąd losowy - model E1")
lines(E1[,2], col='green')
lines(E1[,3], col='blue')
hist(E1, main="Histogram - model E1")

#Error E2
sigmaE2 = matrix(c(1,0,0,0,1,0.9,0,0.9,1), var_size)
E2 = mvrnorm(n=series_length,mu = c(0,0,0), Sigma = sigmaE2)
plot(E2[,1], type='l', col='red', main= "Błąd losowy - model E2")
lines(E2[,2], col='green')
lines(E2[,3], col='blue')
hist(E2, main="Histogram - model E2")

#Error E3
sigma1 = 1
sigma2 = 2
E3.1 = mvrnorm(n=series_length/2,mu = c(0,0,0), Sigma = sigma1*sigma1*diag(3))
E3.2 = mvrnorm(n=series_length/2,mu = c(0,0,0), Sigma = sigma2*sigma2*diag(3))
E3 = rbind(E3.1,E3.2)
plot(E3[,1], type='l', col='red', main= "Błąd losowy - model E3")
lines(E3[,2], col='green')
lines(E3[,3], col='blue')
hist(E3, main="Histogram - model E3")

#Error E4
sigma1 = 1
sigma2 = 3
p = 0.7
E4.1 = mvrnorm(n=series_length,mu = c(0,0,0), Sigma = sigma1*sigma1*diag(3))
E4.2 = mvrnorm(n=series_length,mu = c(0,0,0), Sigma = sigma2*sigma2*diag(3))
E4 = E4.1
s = rbinom(series_length, 1, 0.7)
print(sum(s))
for_i = c(seq(1,series_length))
for(i in for_i)
{
  E.4[i] = s[i]*E4.1[i] + (1-s[i])*E4.2[i]
}
plot(E4[,1], type='l', col='red', main= "Błąd losowy - model E4")
lines(E4[,2], col='green')
lines(E4[,3], col='blue')
hist(E4, main="Histogram - model E4")

#Error E5 - model ARCH
w1 = rnorm(series_length, 0, 1)
w2 = rnorm(series_length, 0, 1)
w3 = rnorm(series_length, 0, 1)
E5 = matrix(rep(0, series_length*var_size), series_length)
for_i = c(seq(3,series_length))
for(i in for_i)
{
  E5[i,1] = w1[i]*sqrt(0.5+(0.1*E5[i-1,1]^2)+(0.4*E5[i-2,1]^2))
  E5[i,2] = w2[i]*sqrt(0.5+(0.1*E5[i-1,2]^2)+(0.4*E5[i-2,2]^2))
  E5[i,3] = w3[i]*sqrt(0.5+(0.1*E5[i-1,3]^2)+(0.4*E5[i-2,3]^2))
  print(E5[i,1])

}
plot(E5[,1], type='l', col='red', main= "Błąd losowy - model E5")
lines(E5[,2], col='green')
lines(E5[,3], col='blue')
hist(E5, main="Histogram - model E5")

par(mfrow=c(3,2))
hist(E1, main="Histogram - błąd losowy - model E1")
hist(E2, main="Histogram - błąd losowy - model E2")
hist(E3, main="Histogram - błąd losowy - model E3")
hist(E4, main="Histogram - błąd losowy - model E4")
hist(E5, main="Histogram - błąd losowy - model E5")

#Funkcja tworzaca model VAR dla okrelonej Macierzy A i Bledu E
create_VAR = function(A, E, series_length, var_size, c_const)
{
  for_i = c(seq(10,series_length))
  YAE = mvrnorm(n=series_length,mu = c(0,0,0), Sigma = diag(3))
  for(i in for_i)
  {
    YAE[i,]=c_const+(A%*%YAE[i-1,])+E[i,] #Model VAR: yt = c+A*yt-1 + et
  }
  return(YAE)
}


bootstrap = function(VAR_model, lagp, N)
{
  #Bootstrap oparty na resztach (z dźwignią)

  #Estymacja za pomocą MNK
  y1 = VAR_model[,1]
  y2 = VAR_model[,2]
  y3 = VAR_model[,3]

  y1_m = lm(y1[2:series_length] ~ y1[1:(series_length-1)] + y3[2:series_length]) #H0
  y2_m = lm(y2[2:series_length] ~ y1[2:series_length] + y2[1:(series_length-1)] + y3[2:series_length])
  y3_m = lm(y3[2:series_length] ~ y1[2:series_length] + y2[2:series_length] + y3[1:(series_length-1)])

  y1_fv = t(as.matrix(y1_m$fitted.values))
  y2_fv = t(as.matrix(y2_m$fitted.values))
  y3_fv = t(as.matrix(y3_m$fitted.values))

  Y1_b = y1_fv
  Y2_b = y2_fv
  Y3_b = y3_fv

  y1_res = y1_m$residuals
  y2_res = y2_m$residuals
  y3_res = y3_m$residuals

  #Odjecie sredniej z reszt
  y1_res = y1_res - mean(y1_res)
  y2_res = y2_res - mean(y2_res)
  y3_res = y3_res - mean(y3_res)

  #Macierz bledow
  y1_res = t(as.matrix(y1_res))
  y2_res = t(as.matrix(y2_res))
  y3_res = t(as.matrix(y3_res))
  Errors = rbind(y1_res,y2_res)
  Errors = rbind(Errors,y3_res)


  for_N = seq(1, N)
  Wyniki_Bootstrap = matrix(0,1,N)

  Errors_resampled = matrix(0,3,series_length)

  for(i in for_N)
  {
    #Losowanie ze zwracaniem + odjecie srednich
    Errors_resampled[1,] = sample(Errors[1,], size=series_length, replace = TRUE)
    Errors_resampled[1,] = Errors_resampled[1,] - mean(Errors_resampled[1,])
    Errors_resampled[2,] = sample(Errors[2,], size=series_length, replace = TRUE)
    Errors_resampled[2,] = Errors_resampled[2,] - mean(Errors_resampled[2,])
    Errors_resampled[3,] = sample(Errors[3,], size=series_length, replace = TRUE)
    Errors_resampled[3,] = Errors_resampled[3,] - mean(Errors_resampled[3,])

    #Wygenerowanie nowych danych z uzyciem oryginalnych danych oraz reszt bootstrapowych
    Y1_b = y1_fv+Errors_resampled[1,1:(series_length-1)]
    Y2_b = y2_fv+Errors_resampled[2,1:(series_length-1)]
    Y3_b = y3_fv+Errors_resampled[3,1:(series_length-1)]
    VAR_model_b = rbind(Y1_b,Y2_b)
    VAR_model_b = rbind(VAR_model_b,Y3_b)
    VAR_model_b = t(VAR_model_b)

    #Estymacja modelu VAR - jako stopień opóżnienia: parametr opoznienia + maksymalny stopień integracji (Zawsze zakładamy 1)
    d = 1 #Order of integration - 1
    V1 = VAR(VAR_model_b[101:series_length-1,],p=lagp+d,type="const") #VAR method is using ols (MNK)
    summary(V1)

    #Macierz beta
    Beta=rbind(as.matrix(V1$varresult$y1$coefficients),as.matrix(V1$varresult$y2$coefficients))
    Beta=as.matrix(rbind(Beta,as.matrix(V1$varresult$y3$coefficients)))

    #macierz wariancji-kowariancji
    SigmaU = summary(V1)$covres

    #Macierz Y
    Y = t(VAR_model_b)

    #Macierz Z
    Z = matrix(0,1+var_size*(lagp+d),series_length-99)
    for(k in 1:(series_length-99))
    {
      dd=matrix(1,1,1)
      for (l in 0:lagp)
      {
        Y_pom=matrix(Y[,(99-l+k-1)],3,1)
        dd=rbind(dd,Y_pom)
      }
      Z[,k]=dd
    }

    #Macierz C (Macierz 0 lub 1 o wymiarze p x (1+n(p+d)))
    C=matrix(0,1,var_size*(1+var_size*(lagp+d)))
    #Ustalenie hipotezy 0: y2 nie jest przyczyną w sensie Grangera dla y1
    C[4]=1

    #obliczenie statystyki TY (testu Toda-Yamamoto)
    TY_stat<-t(C%*%Beta)%*%solve(C%*%(kronecker(solve(Z%*%t(Z)),SigmaU))%*%t(C))%*%(C%*%Beta)
    TY2<-TY_stat[1,1]

    Wyniki_Bootstrap[1,i]<-TY2
  }


  return(Wyniki_Bootstrap)
}

#_____________GŁÓWNA PĘTŁA__________________________________________________________________________

#Funkcja do obliczania ile razy test odrzuci H0
rejection = function(lagp, number_of_repetitions, N, A, E, series_length, var_size, c_const)
{

  #Statystyka Chi^2 dla konkretnych poziomów istotności i jednego stopnia swobody:
  chi0.01 = 6.635
  chi0.05 = 3.841
  chi0.1 = 2.706

  #Poziomy istotności dla bootstrapu
  poziom_b0.1 = 0.9
  poziom_b0.05 = 0.95
  poziom_b0.01 = 0.99

  number_of_rejections0.01 = 0 #Zmienna kontrolujaca ilosc odrzucen H0 przez test na poziomie istotnosci 0.01
  number_of_rejections0.05 = 0 #Zmienna kontrolujaca ilosc odrzucen H0 przez test na poziomie istotnosci 0.05
  number_of_rejections0.1 = 0 #Zmienna kontrolujaca ilosc odrzucen H0 przez test na poziomie istotnosci 0.1
  number_of_rejections0.01_b = 0 #Zmienna kontrolujaca ilosc odrzucen H0 przez test na poziomie istotnosci 0.01 przy zastosowaniu bootstrapu
  number_of_rejections0.05_b = 0 #Zmienna kontrolujaca ilosc odrzucen H0 przez test na poziomie istotnosci 0.05 przy zastosowaniu bootstrapu
  number_of_rejections0.1_b = 0 #Zmienna kontrolujaca ilosc odrzucen H0 przez test na poziomie istotnosci 0.1 przy zastosowaniu bootstrapu


  for(j in 1:number_of_repetitions)
  {
    #Tworzenie modelu VAR
    if(E == 'E1')
    {
      #Error E1
      E1 = mvrnorm(n=series_length,mu = c(0,0,0), Sigma = diag(3))
      Er = E1
    }
    else if(E == 'E2')
    {
      #Error E2
      sigmaE2 = matrix(c(1,0,0,0,1,0.9,0,0.9,1), var_size)
      E2 = mvrnorm(n=series_length,mu = c(0,0,0), Sigma = sigmaE2)
      Er = E2
    }
    else if(E == 'E3')
    {
      #Error E3
      sigma1 = 1
      sigma2 = 2
      E3.1 = mvrnorm(n=series_length/2,mu = c(0,0,0), Sigma = sigma1*sigma1*diag(3))
      E3.2 = mvrnorm(n=series_length/2,mu = c(0,0,0), Sigma = sigma2*sigma2*diag(3))
      E3 = rbind(E3.1,E3.2)
      Er = E3
    }
    else if(E == 'E4')
    {
      #Error E4
      sigma1 = 1
      sigma2 = 3
      p = 0.7
      E4.1 = mvrnorm(n=series_length,mu = c(0,0,0), Sigma = sigma1*sigma1*diag(3))
      E4.2 = mvrnorm(n=series_length,mu = c(0,0,0), Sigma = sigma2*sigma2*diag(3))
      E4 = E4.1
      s = rbinom(series_length, 1, 0.7)
      for_i = c(seq(1,series_length))
      for(i in for_i)
      {
        E4[i] = s[i]*E4.1[i] + (1-s[i])*E4.2[i]
      }
      Er = E4
    }
    else if(E == 'E5')
    {
      #Error E5 - model ARCH
      w1 = rnorm(series_length, 0, 1)
      w2 = rnorm(series_length, 0, 1)
      w3 = rnorm(series_length, 0, 1)
      E5 = matrix(rep(0, series_length*var_size), series_length)
      for_i = c(seq(3,series_length))
      for(i in for_i)
      {
        E5[i,1] = w1[i]*sqrt(0.5+(0.1*E5[i-1,1]^2)+(0.4*E5[i-2,1]^2))
        E5[i,2] = w2[i]*sqrt(0.5+(0.1*E5[i-1,2]^2)+(0.4*E5[i-2,2]^2))
        E5[i,3] = w3[i]*sqrt(0.5+(0.1*E5[i-1,3]^2)+(0.4*E5[i-2,3]^2))
      }
      Er = E5
    }
    VAR_model = create_VAR(A, Er, series_length, var_size, c_const)

    #Estymacja modelu VAR - jako stopień opóżnienia: parametr opoznienia + maksymalny stopień integracji (Zawsze zakładamy 1)
    d = 1 #Order of integration - 1
    V1 = VAR(VAR_model[101:series_length,],p=lagp+d,type="const") #VAR method is using ols (MNK)
    summary(V1)

    #Macierz beta
    Beta=rbind(as.matrix(V1$varresult$y1$coefficients),as.matrix(V1$varresult$y2$coefficients))
    Beta=as.matrix(rbind(Beta,as.matrix(V1$varresult$y3$coefficients)))

    #macierz wariancji-kowariancji
    SigmaU = summary(V1)$covres

    #Macierz Y
    Y = t(VAR_model)

    #Macierz Z
    Z = matrix(0,1+var_size*(lagp+d),series_length-100)
    for(k in 1:(series_length-100))
    {
      dd=matrix(1,1,1)
      for (l in 0:lagp)
      {
        Y_pom=matrix(Y[,(100-l+k-1)],3,1)
        dd=rbind(dd,Y_pom)
      }
      Z[,k]=dd
    }

    #Macierz C (Macierz 0 lub 1 o wymiarze p x (1+n(p+d)))
    C=matrix(0,1,var_size*(1+var_size*(lagp+d)))
    #Ustalenie hipotezy 0: y2 nie jest przyczyną w sensie Grangera dla y1
    C[4]=1

    #obliczenie statystyki TY (testu Toda-Yamamoto)
    TY_stat<-t(C%*%Beta)%*%solve(C%*%(kronecker(solve(Z%*%t(Z)),SigmaU))%*%t(C))%*%(C%*%Beta)
    TY<-TY_stat[1,1]

    # result = wald.test(b=coef(V1$varresult$y1), Sigma=vcov(V1$varresult$y1), Terms = 2)
    # TY=result$result$chi2[1]

    if(TY > chi0.01){
      number_of_rejections0.01 = number_of_rejections0.01+1
    }

    if(TY > chi0.05)
    {
      number_of_rejections0.05 = number_of_rejections0.05+1
    }

    if(TY > chi0.1)
    {
      number_of_rejections0.1 = number_of_rejections0.1+1
    }


    #____________________________________________________________________________________

    Wyniki_Bootstrap = bootstrap(VAR_model, lagp, N)
    Wyniki_Bootstrap = t(as.matrix(sort(Wyniki_Bootstrap)))
    Wyniki_Bootstrap2 = as.matrix(unique(Wyniki_Bootstrap))
    
    dystrybuanta = matrix(1,1,N) #Dystrybuanta dla uzyskania kwantylu

    ld=ncol(dystrybuanta)
    dystrybuanta=t(as.matrix(dystrybuanta[1,2:ld]))
    ld=ld-1

    for (id in 1:(ld-1)) {
      dystrybuanta[1,id+1]=dystrybuanta[1,id]+1
    }

    dystrybuanta=(1/N)*dystrybuanta

    #Ustalenie poziomu dla różnych poziomów istotności

    przedzial0.1=0

    for (i0.1 in 1:(ld-1)) {
      if (dystrybuanta[i0.1]<=poziom_b0.1) {
        if (dystrybuanta[i0.1+1]>poziom_b0.1) {
          przedzial0.1=i0.1
        }
      }
    }
    
    przedzial0.05=0
    
    for (i0.05 in 1:(ld-1)) {
      if (dystrybuanta[i0.05]<=poziom_b0.05) {
        if (dystrybuanta[i0.05+1]>poziom_b0.05) {
          przedzial0.05=i0.05
        }
      }
    }
    
    przedzial0.01=0
    
    for (i0.01 in 1:(ld-1)) {
      if (dystrybuanta[i0.01]<=poziom_b0.01) {
        if (dystrybuanta[i0.01+1]>poziom_b0.01) {
          przedzial0.01=i0.01
        }
      }
    }
    if(przedzial0.01 == 0)
    {
      przedzial0.01 = ld
    }

    
    waga0.1=0
    if (poziom_b0.1-dystrybuanta[przedzial0.1]<0) {
      waga0.1=(dystrybuanta[przedzial0.1+1]-poziom_b0.1)/(poziom_b0.1-dystrybuanta[przedzial0.1])
    }
    kwantyl0.1=Wyniki_Bootstrap2[przedzial0.1]+waga0.1*Wyniki_Bootstrap2[przedzial0.1+1]
    if (TY>kwantyl0.1) {
      number_of_rejections0.1_b=number_of_rejections0.1_b+1
    }


    waga0.05=0
    if (poziom_b0.05-dystrybuanta[przedzial0.05]<0) {
      waga0.05=(dystrybuanta[przedzial0.05+1]-poziom_b0.05)/(poziom_b0.05-dystrybuanta[przedzial0.05])
    }
    kwantyl0.05=Wyniki_Bootstrap2[przedzial0.05]+waga0.05*Wyniki_Bootstrap2[przedzial0.05+1]
    if (TY>kwantyl0.05) {
      number_of_rejections0.05_b=number_of_rejections0.05_b+1
    }


    waga0.01=0
    if (poziom_b0.01-dystrybuanta[przedzial0.01]<0) {
      waga0.01=(dystrybuanta[przedzial0.01+1]-poziom_b0.01)/(poziom_b0.01-dystrybuanta[przedzial0.01])
    }
    kwantyl0.01=Wyniki_Bootstrap2[przedzial0.01]+waga0.01*Wyniki_Bootstrap2[przedzial0.01+1]
    if (TY>kwantyl0.01) {
      number_of_rejections0.01_b=number_of_rejections0.01_b+1
    }

    #____________________________________________________________________________________

  }
  wektor_odrzucen = c(number_of_rejections0.01/number_of_repetitions,number_of_rejections0.05/number_of_repetitions,number_of_rejections0.1/number_of_repetitions, number_of_rejections0.01_b/number_of_repetitions,number_of_rejections0.05_b/number_of_repetitions,number_of_rejections0.1_b/number_of_repetitions)
  return(wektor_odrzucen)
  
}

# #Dane A1 + E1
# YA1E1 = create_VAR(A1, E1, series_length, var_size, c_const)
# #Dane A2 + E1
# YA2E1 = create_VAR(A2, E1, series_length, var_size, c_const)
# #Dane A3 + E1
# YA3E1 = create_VAR(A3, E1, series_length, var_size, c_const)
#
# par(mfrow = c(1, 1))
# plot(YA1E1[,1], type='l', col='red', ylim=range(-15:15), main="Model A1 - blad E1")
# lines(YA1E1[,2], col='green')
# lines(YA1E1[,3], col='blue')
#
# plot(YA2E1[,1], type='l', col='red', ylim=range(-15:15), main="Model A2 - blad E1")
# lines(YA2E1[,2], col='green')
# lines(YA2E1[,3], col='blue')
#
# plot(YA3E1[,1], type='l', col='red', ylim=range(-15:15), main="Model A3 - blad E1")
# lines(YA3E1[,2], col='green')
# lines(YA3E1[,3], col='blue')

#HYPOTESIS TO TEST: y2 do not Granger-cause y1

lagp = 1 #Wartosc parametru p = 1 (Dobrze okreslony parametr opoznien) lub p = 2 (Zle okreslony parametr opoznien)
number_of_repetitions = 100 #Ilosc powtorzen testu
N = 1000 #Ilość realizacji bootstrapowych


print("Percent of rejections for: A1 + E1")
wod = rejection(lagp, number_of_repetitions, N, A1, 'E1', series_length, var_size, c_const)
print("1%: " + wod[1] + ", 5%: " + wod[2] + ", 10%: " + wod[3] + ", (B)1%: " + wod[4]+ ", (B)5%: " + wod[5] + ", (B)10%: " +wod[6])
