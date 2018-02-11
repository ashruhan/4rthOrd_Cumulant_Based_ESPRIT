mm=[1 1
    2 2
    3 3
    1 2
    1 3
    2 3];

R10 = p1*p1.'/Window;
R20 = p2*p2.'/Window;

 for m=1:6
    
            for n=1:6
        
                R11(m,n)= R11(m,n) - R10(mm(m,1),mm(m,2))*R10(mm(n,1),mm(n,2));
                R11(m,n)= R11(m,n) - R1(mm(m,1),mm(n,1))*R1(mm(m,2),mm(n,2));
                R11(m,n)= R11(m,n) - R1(mm(m,1),mm(n,2))*R1(mm(m,2),mm(n,1));
                
                R22(m,n)= R22(m,n) - R10(mm(m,1),mm(m,2))*R20(mm(n,1),mm(n,2));
                R22(m,n)= R22(m,n) - R2(mm(m,1),mm(n,1))*R2(mm(m,2),mm(n,2));
                R22(m,n)= R22(m,n) - R2(mm(m,1),mm(n,2))*R2(mm(m,2),mm(n,1));

            end

 end