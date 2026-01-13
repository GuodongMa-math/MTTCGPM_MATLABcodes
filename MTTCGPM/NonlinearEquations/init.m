% init 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [nprob,n,m,x0]=initf(nprob)
% This function sets n,m, and the standard starting    
% point based on the nprob and returns it to initpt     
% function.                                                                                                    
% Created on 10/30/94 by Madhu Lamba                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nprob,x0,x0_old,x0_old1,n] = init(NO)

%global FIRSTIME

switch NO
    %% nprob='problem1'; cuter
      case 1
        nprob='problem1';
        n=10000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
        
        
    %% nprob='problem1'; cuter
      case 2
        nprob='problem1';
        n=10000;
        x0=0.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
       
    %% nprob='problem1'; cuter
      case 3
        nprob='problem1';
        n=10000;
        x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 
    %% nprob ='problem1'; cuter
      case 4
        nprob='problem1';
        n=10000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 
   
    %% nprob='problem1'; cuter
      case 5
        nprob='problem1';
        n=10000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1)./0.2; 
    %% nprob='problem1'; cuter
      case 6
        nprob='problem1';
        n=10000;
         x0=2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 
   case 7
        nprob='problem1';
        n=10000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);      
             
    %% nprob='problem1'; cuter
    case 8
        nprob='problem1';
        n=50000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
        %x0=ones(n,1);
    
    %% nprob='problem1'; cuter
      case 9
        nprob='problem1';
        n=50000;
        x0=0.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 10
        nprob='problem1';
        n=50000;
       x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
       
    %% nprob='problem1'; cuter
      case 11
        nprob='problem1';
        n=50000;
      x0=1.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 12
        nprob='problem1';
        n=50000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 13
        nprob='problem1';
        n=50000;
       x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
      
    %% nprob='problem1'; cuter
      case 14
        nprob='problem1';
        n=50000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);  

        case 15
        nprob='problem1';
        n=100000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
           

    %% nprob='problem1'; cuter
      case 16
        nprob='problem1';
        n=100000;
         x0=0.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
  
    %% nprob='problem1'; cuter
      case 17
        nprob='problem1';
        n=100000;
          x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
       
    %% nprob='problem1'; cuter
      case 18
        nprob='problem1';
        n=100000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
       
    %% nprob='problem1'; cuter
      case 19
        nprob='problem1';
        n=100000;
         x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
  
    %% nprob='problem1'; cuter
      case 20
        nprob='problem1';
        n=100000;
         x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

      case 21
        nprob='problem1';
        n=100000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);        
        
    %% nprob='problem1'; cuter
    case 22
        nprob='problem2';
        n=10000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
               
    %% nprob='problem1'; cuter
      case 23
        nprob='problem2';
        n=10000;
        x0=0.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 24
        nprob='problem2';
        n=10000;
        x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 
        
    %% nprob='problem1'; cuter
      case 25
        nprob='problem2';
        n=10000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 
    
    %% nprob='problem1'; cuter
      case 26
        nprob='problem2';
        n=10000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 
        
    %% nprob='problem1'; cuter
      case 27
        nprob='problem2';
        n=10000;
         x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 

      case 28
        nprob='problem2';
        n=10000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);        

    %% nprob='problem1'; cuter
    case 29
        nprob='problem2';
        n=50000;
        x0=0.1*ones(n,1);
       x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
     
    %% nprob='problem1'; cuter
      case 30
        nprob='problem2';
        n=50000;
        x0=0.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 31
        nprob='problem2';
        n=50000;
       x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 32
        nprob='problem2';
        n=50000;
      x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
  
    %% nprob='problem1'; cuter
      case 33
        nprob='problem2';
        n=50000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 34
        nprob='problem2';
        n=50000;
       x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 35
        nprob='problem2';
        n=50000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

      case 36
        nprob='problem2';
        n=100000;
         x0=0.1*ones(n,1);
       x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);        
    %% nprob='problem1'; cuter
      case 37
        nprob='problem2';
        n=100000;
         x0=0.2*ones(n,1);
        x0_old=0.1*ones(n,1);
       x0_old1=0.2*ones(n,1);
        %x0=ones(n,1);
    
    %% nprob='problem1'; cuter
      case 38
        nprob='problem2';
        n=100000;
          x0=0.5*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 39
        nprob='problem2';
        n=100000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 40
        nprob='problem2';
        n=100000;
         x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
   
    %% nprob='problem1'; cuter
      case 41
        nprob='problem2';
        n=100000;
         x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        case 42
        nprob='problem2';
        n=100000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);      
        
    %% nprob='problem1'; cuter   
       case 64
        nprob='problem3';
        n=10000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
        
        
    %% nprob='problem1'; cuter
      case 65
        nprob='problem3';
        n=10000;
        x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
       
    %% nprob='problem1'; cuter
      case 66
        nprob='problem3';
        n=10000;
        x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
    %% nprob='problem1'; cuter
      case 67
        nprob='problem3';
        n=10000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 
    
    %% nprob='problem1'; cuter
      case 68
        nprob='problem3';
        n=10000;
       x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
    %% nprob='problem1'; cuter
      case 69
        nprob='problem3';
        n=10000;
         x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

        case 70
        nprob='problem3';
        n=10000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);       
    
    %% nprob='problem1'; cuter
    case 71
        nprob='problem3';
        n=50000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
      
    %% nprob='problem1'; cuter
      case 72
        nprob='problem3';
        n=50000;
        x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 73
        nprob='problem3';
        n=50000;
       x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 74
        nprob='problem3';
        n=50000;
      x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
   
    %% nprob='problem1'; cuter
      case 75
        nprob='problem3';
        n=50000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 76
        nprob='problem3';
        n=50000;
       x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 77
        nprob='problem3';
        n=50000;
        x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 

        case 78
        nprob='problem3';
        n=100000;
         x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);        
    %% nprob='problem1'; cuter
      case 79
        nprob='problem3';
        n=100000;
         x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
      
    %% nprob='problem1'; cuter
      case 80
        nprob='problem3';
        n=100000;
          x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 81
        nprob='problem3';
        n=100000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 82
        nprob='problem3';
        n=100000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 83
        nprob='problem3';
        n=100000;
         x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
    case 84
        nprob='problem3';
        n=100000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);        
    case 85
        nprob='problem4';
        n=10000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
        
        
    %% nprob='problem1'; cuter
      case 86
        nprob='problem4';
        n=10000;
        x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 87
        nprob='problem4';
        n=10000;
        x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
    %% nprob='problem1'; cuter
      case 88
        nprob='problem4';
        n=10000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
    %% nprob='problem1'; cuter
      case 89
        nprob='problem4';
        n=10000;
       x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
    %% nprob='problem1'; cuter
      case 90
        nprob='problem4';
        n=10000;
        x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

        case 91
        nprob='problem4';
        n=10000;
          x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);         
 
    %% nprob='problem1'; cuter
    case 92
        nprob='problem4';
        n=50000;
       x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
        %x0=ones(n,1);
    
    %% nprob='problem1'; cuter
      case 93
        nprob='problem4';
        n=50000;
       x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 94
        nprob='problem4';
        n=50000;
       x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 95
        nprob='problem4';
        n=50000;
      x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 96
        nprob='problem4';
        n=50000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 97
        nprob='problem4';
        n=50000;
       x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 98
        nprob='problem4';
        n=50000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);  

        case 99
        nprob='problem4';
        n=100000;
         x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
             
    %% nprob='problem1'; cuter
      case 100
        nprob='problem4';
        n=100000;
         x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
     
    %% nprob='problem1'; cuter
      case 101
        nprob='problem4';
        n=100000;
         x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 102
        nprob='problem4';
        n=100000;
         x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 103
        nprob='problem4';
        n=100000;
       x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 104
        nprob='problem4';
        n=100000;
        x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
     case 105
        nprob='problem4';
        n=100000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);       
          
    %% nprob='problem1'; cuter       
         case 127
        nprob='problem5';
        n=10000;
      x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
        
        
    %% nprob='problem1'; cuter
      case 128
        nprob='problem5';
        n=10000;
        x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
       
    %% nprob='problem1'; cuter
      case 129
        nprob='problem5';
        n=10000;
       x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
    %% nprob='problem1'; cuter
      case 130
        nprob='problem5';
        n=10000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 131
        nprob='problem5';
        n=10000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
    %% nprob='problem1'; cuter
      case 132
        nprob='problem5';
        n=10000;
         x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

      case 133
        nprob='problem5';
        n=10000;
          x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);         
  
    %% nprob='problem1'; cuter
    case 134
        nprob='problem5';
        n=50000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
        %x0=ones(n,1);
    
    %% nprob='problem1'; cuter
      case 135
        nprob='problem5';
        n=50000;
        x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 136
        nprob='problem5';
        n=50000;
       x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 137
        nprob='problem5';
        n=50000;
      x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 138
        nprob='problem5';
        n=50000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 139
        nprob='problem5';
        n=50000;
       x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 140
        nprob='problem5';
        n=50000;
          x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);   

    %% nprob='problem1'; cuter
      case 141
        nprob='problem5';
        n=100000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
 
    %% nprob='problem1'; cuter
      case 142
        nprob='problem5';
        n=100000;
        x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
  
    %% nprob='problem1'; cuter
      case 143
        nprob='problem5';
        n=100000;
         x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
  
    %% nprob='problem1'; cuter
      case 144
        nprob='problem5';
        n=100000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
 
    %% nprob='problem1'; cuter
      case 145
        nprob='problem5';
        n=100000;
         x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
         
      case 146
        nprob='problem5';
        n=100000;
         x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
                 
      case 147
        nprob='problem5';
        n=100000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 
       
        
        
    %% nprob='problem1'; cuter
      case 148
        nprob='problem6';
        n=10000;
       x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);

    %% nprob='problem1'; cuter
      case 149
        nprob='problem6';
        n=10000;
        x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 150
        nprob='problem6';
        n=10000;
        x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 

    %% nprob='problem1'; cuter
      case 151
        nprob='problem6';
        n=10000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 
       
    %% nprob='problem1'; cuter
      case 152
        nprob='problem6';
        n=10000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
        case 153
        nprob='problem6';
        n=10000;
         x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        case 154
        nprob='problem6';
        n=10000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);        
 
    %% nprob='problem1'; cuter
    case 155
        nprob='problem6';
        n=50000;
       x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
        %x0=ones(n,1);
    
    %% nprob='problem1'; cuter
      case 156
        nprob='problem6';
        n=50000;
        x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 157
        nprob='problem6';
        n=50000;
        x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 158
        nprob='problem6';
        n=50000;
       x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 159
        nprob='problem6';
        n=50000;
       x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 160
        nprob='problem6';
        n=50000;
       x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 161
        nprob='problem6';
        n=50000;
       x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 

    %% nprob='problem1'; cuter
      case 162
        nprob='problem6';
        n=100000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
       
  
    
    %% nprob='problem1'; cuter
      case 163
        nprob='problem6';
        n=100000;
         x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
          
    %% nprob='problem1'; cuter
      case 164
        nprob='problem6';
        n=100000;
        x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        

    %% nprob='problem1'; cuter
      case 165
        nprob='problem6';
        n=100000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 166
        nprob='problem6';
        n=100000;
         x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
        case 167
        nprob='problem6';
        n=100000;
         x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        case 168
        nprob='problem6';
        n=100000;
          x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);       
         case 169
        nprob='problem7';
        n=10000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
        
        
    %% nprob='problem1'; cuter
      case 170
        nprob='problem7';
        n=10000;
        x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 171
        nprob='problem7';
        n=10000;
      x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 172
        nprob='problem7';
        n=10000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
  
    %% nprob='problem1'; cuter
      case 173
        nprob='problem7';
        n=10000;
       x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
    %% nprob='problem1'; cuter
      case 174
        nprob='problem7';
        n=10000;
         x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
    

    case 175
        nprob='problem7';
        n=10000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);         
    %% nprob='problem1'; cuter
    case 176
        nprob='problem7';
        n=50000;
       x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
        %x0=ones(n,1);
    
    %% nprob='problem1'; cuter
      case 177
        nprob='problem7';
        n=50000;
       x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 178
        nprob='problem7';
        n=50000;
      x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 179
        nprob='problem7';
        n=50000;
      x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 180
        nprob='problem7';
        n=50000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 181
        nprob='problem7';
        n=50000;
       x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 182
        nprob='problem7';
        n=50000;
        x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
         
    case 183
        nprob='problem7';
        n=100000;
         x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);        
    %% nprob='problem1'; cuter
      case 184
        nprob='problem7';
        n=100000;
       x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 185
        nprob='problem7';
        n=100000;
          x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 186
        nprob='problem7';
        n=100000;
         x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 187
        nprob='problem7';
        n=100000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 188
        nprob='problem7';
        n=100000;
         x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
     case 189
        nprob='problem7';
        n=100000;
          x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);      

        case 209
        nprob='problem8';
        n=10000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);

        %% nprob='problem1'; cuter
      case 210
        nprob='problem8';
        n=10000;
        x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
       
    %% nprob='problem1'; cuter
      case 211
        nprob='problem8';
        n=10000;
        x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 212
        nprob='problem8';
        n=10000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
      
    %% nprob='problem1'; cuter
      case 213
        nprob='problem8';
        n=10000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
  
    case 214
        nprob='problem8';
        n=10000;
        x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
         case 215
        nprob='problem8';
        n=10000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);       
    %% nprob='problem1'; cuter
    case 216
        nprob='problem8';
        n=50000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);
    
    %% nprob='problem1'; cuter
      case 217
        nprob='problem8';
        n=50000;
        x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 218
        nprob='problem8';
        n=50000;
        x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 219
        nprob='problem8';
        n=50000;
      x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 220
        nprob='problem8';
        n=50000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 221
        nprob='problem8';
        n=50000;
      x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
    %% nprob='problem1'; cuter
      case 222
        nprob='problem8';
        n=50000;
        x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);  
           
    %% nprob='problem1'; cuter
      case 223
        nprob='problem8';
        n=100000;
        x0=0.1*ones(n,1);
        x0_old=0.2*ones(n,1);
        x0_old1=0.1*ones(n,1);      
    
    %% nprob='problem1'; cuter
      case 224
        nprob='problem8';
        n=100000;
        x0=0.2*ones(n,1);
       x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 225
        nprob='problem8';
        n=100000;
        x0=0.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
  
    %% nprob='problem1'; cuter
      case 226
        nprob='problem8';
        n=100000;
        x0=1.2*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);

    %% nprob='problem1'; cuter
      case 227
        nprob='problem8';
        n=100000;
        x0=1.5*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
        
        case 228
        nprob='problem8';
        n=100000;
         x0=2.0*ones(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1);
                
     case 229
        nprob='problem8';
        n=100000;
         x0=rand(n,1);
        x0_old=0.1*ones(n,1);
        x0_old1=0.2*ones(n,1); 

end