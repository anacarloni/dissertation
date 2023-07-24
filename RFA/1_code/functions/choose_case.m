% Setup do Caso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch caso
    case '01'
        nsample = 2;
        novlp = 0;
        window = 'hanning';
    case '02'
        nsample = 1;
        novlp = 0;
        window = 'rectwin';
    case '03'
        nsample = 3;
        novlp = 0;
        window = 'hanning';
    case '04'
        nsample = 2;
        novlp = 1;
        window = 'rectwin';
    case '05'
        nsample = 1;
        novlp = 0;
        window = 'rectwin';
    case '06'
        nsample = 3;
        novlp = 0;
        window = 'hanning';
    case '07'
        nsample = 2;
        novlp = 1;
        window = 'hanning';
    case '08'
        nsample = 1;
        novlp = 0;
        window = 'rectwin';
    case '09'
        nsample = 4;
        novlp = 0;
        window = 'hanning';
    case '10'
        nsample = 3;
        novlp = 2;
        window = 'hanning';
    case '11'
        nsample = 2;
        novlp = 1;
        window = 'hanning';
    case '12'
        nsample = 2;
        novlp = 0;
        window = 'hanning';
    case '13'
        nsample = 1;
        novlp = 0;
        window = 'rectwin';
    case '14'
        nsample = 4;
        novlp = 0;
        window = 'hanning';
    case '15'
        nsample = 3;
        novlp = 2;
        window = 'hanning';
    case '16'
        nsample = 2;
        novlp = 1;
        window = 'hanning';
    case '17'
        nsample = 2;
        novlp = 0;
        window = 'hanning';
    case '18'
        nsample = 1;
        novlp = 0;
        window = 'rectwin';
    case '19'
        nsample = 4;
        novlp = 0;
    case '20'
        nsample = 3;
        novlp = 2;
    case '21'
        nsample = 2;
        novlp = 1;
    case '22'
        nsample = 2;
        novlp = 0;
    case '23'
        nsample = 1;
        novlp = 0;
    case '24'
        nsample = 4;
        novlp = 0;
    case '25'
        nsample = 3;
        novlp = 2;
    case '26'
        nsample = 2;
        novlp = 1;
    case '27'
        nsample = 2;
        novlp = 0;
    case '28'
        nsample = 1;
        novlp = 0;
end