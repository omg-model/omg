function [ gchem_fcns ] = geochemistry_functions ( )
    
    gchem_fcns.calc_constants_CO2SYS=@calc_constants_CO2SYS;
    gchem_fcns.calc_constants_OMG=@calc_constants_OMG;
    gchem_fcns.interpolate_constants=@interpolate_constants;
    gchem_fcns.calc_carbonate_constants=@calc_carbonate_constants;
    gchem_fcns.calc_constants=@calc_constants;
    gchem_fcns.solve_carbonate_system=@solve_carbonate_system;

end

%%----- SUBROUTINES -----%%
%%
function [constants] = calc_constants_CO2SYS(constants,I,T,P,S,pHScale,K1K2,KSO4)

    % Carbonate chemistry constants, based on CO2SYS Constants function, modified with globals removed
    % Results are returned in struct C
    %
    %
    % Required input:
    % TempC  - temperature deg C
    % Pdbar  - pressure, dbar 
    % S      - salinity 
    % Optional input:
    % pHScale       - default 3 (free - required for CalculateTAfromTCpHfree)
    % K1K2CONSTANTS - default 10 
    % KSO4CONSTANTS - default 1          
    %
    % Set pH scale
    %  1 = Total scale
    %  2 = Seawater scale
    %  3 = Free scale   [required for copse_CO2SYS]
    %  4 = NBS scale
    %
    %   K1K2CONSTANTS     : scalar or vector of size n (***)
    %  (***) Each element must be an integer,
    %        indicating the K1 K2 dissociation constants that are to be used:
    %   1 = Roy, 1993											T:    0-45  S:  5-45. Total scale. Artificial seawater.
    %   2 = Goyet & Poisson										T:   -1-40  S: 10-50. Seaw. scale. Artificial seawater.
    %   3 = HANSSON              refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
    %   4 = MEHRBACH             refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
    %   5 = HANSSON and MEHRBACH refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
    %   6 = GEOSECS (i.e., original Mehrbach)					T:    2-35  S: 19-43. NBS scale.   Real seawater.
    %   7 = Peng	(i.e., originam Mehrbach but without XXX)	T:    2-35  S: 19-43. NBS scale.   Real seawater.
    %   8 = Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)	T:    0-50  S:     0.
    %   9 = Cai and Wang, 1998									T:    2-35  S:  0-49. NBS scale.   Real and artificial seawater.
    %  10 = Lueker et al, 2000									T:    2-35  S: 19-43. Total scale. Real seawater.
    %  11 = Mojica Prieto and Millero, 2002.					T:    0-45  S:  5-42. Seaw. scale. Real seawater
    %  12 = Millero et al, 2002									T: -1.6-35  S: 34-37. Seaw. scale. Field measurements.
    %  13 = Millero et al, 2006									T:    0-50  S:  1-50. Seaw. scale. Real seawater.
    %  14 = Millero et al, 2010									T:    0-50  S:  1-50. Seaw. scale. Real seawater.
    %
    %  (****) Each element must be an integer that
    %         indicates the KSO4 dissociation constants that are to be used,
    %         in combination with the formulation of the borate-to-salinity ratio to be used.
    %         Having both these choices in a single argument is somewhat awkward,
    %         but it maintains syntax compatibility with the previous version.
    %  1 = KSO4 of Dickson & TB of Uppstrom 1979  (PREFERRED)
    %  2 = KSO4 of Khoo    & TB of Uppstrom 1979
    %  3 = KSO4 of Dickson & TB of Lee 2010
    %  4 = KSO4 of Khoo    & TB of Lee 2010

    nd = size(T);

    TempC               = T;      % vector
    Sal                 = S;        % vector
    sqrSal              = sqrt(Sal);    % vector
    Pbar                = P./10;       % vector
    WhichKs             = K1K2;    % scalar
    WhoseKSO4           = KSO4;  % scalar
    pHScale             = pHScale;    % scalar
    RGasConstant     	= 83.1451;  % ml bar-1 K-1 mol-1, DOEv2
    %RGasConstant     	= 83.14472; % ml bar-1 K-1 mol-1, DOEv3
    %C.RGasConstant      = RGasConstant;
    TempK               = TempC + 273.15;
    RT                  = RGasConstant.*TempK;
    logTempK            = log(TempK);

    volkg               = (1 - 0.001005.*Sal);  % l / kg  volume of 1 kg sw at 1 atm
    
    % SUB Constants, version 04.01, 10-13-97, written by Ernie Lewis.
    % Inputs: pHScale%, WhichKs%, WhoseKSO4%, Sali, TempCi, Pdbar
    % Outputs: K0, K(), T(), fH, FugFac, VPFac
    % This finds the Constants of the CO2 system in seawater or freshwater,
    % corrects them for pressure, and reports them on the chosen pH scale.
    % The process is as follows: the Constants (except KS, KF which stay on the
    % free scale - these are only corrected for pressure) are
    %       1) evaluated as they are given in the literature
    %       2) converted to the SWS scale in mol/kg-SW or to the NBS scale
    %       3) corrected for pressure
    %       4) converted to the SWS pH scale in mol/kg-SW
    %       5) converted to the chosen pH scale
    %
    %       PROGRAMMER'S NOTE: all logs are log base e
    %       PROGRAMMER'S NOTE: all Constants are converted to the pH scale
    %               pHScale% (the chosen one) in units of mol/kg-SW
    %               except KS and KF are on the free scale
    %               and KW is in units of (mol/kg-SW)^2


    % CalculateTB - Total Borate:
    switch WhichKs
        case 8
            % Pure water case.
            TB = zeros(nd);
        case {6,7}       
            TB = 0.0004106.*Sal./35; % in mol/kg-SW
            % this is .00001173.*Sali
            % this is about 1% lower than Uppstrom's value
            % Culkin, F., in Chemical Oceanography,
            % ed. Riley and Skirrow, 1965:
            % GEOSECS references this, but this value is not explicitly
            % given here
        otherwise
            switch WhoseKSO4
                case {1,2} % If user opted for Uppstrom's values:
                    % Uppstrom, L., Deep-Sea Research 21:161-162, 1974:
                    % this is .000416.*Sali./35. = .0000119.*Sali
                    % TB(FF) = (0.000232./10.811).*(Sal(FF)./1.80655); % in mol/kg-SW
                    TB =  0.0004157.*Sal./35; % in mol/kg-SW
                case {3,4} % If user opted for the new Lee values:

                    % Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.
                    % Geochimica Et Cosmochimica Acta 74 (6): 1801?1811.
                    TB =  0.0004326.*Sal./35; % in mol/kg-SW
                otherwise
                    error('unrecognized WhoseKSO4 %g',WhoseKSO4);
            end
    end

    % CalculateTF;
    % Riley, J. P., Deep-Sea Research 12:219-220, 1965:
    % this is .000068.*Sali./35. = .00000195.*Sali
    TF = (0.000067./18.998).*(Sal./1.80655); % in mol/kg-SW

    % CalculateTS ;
    % Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
    % this is .02824.*Sali./35. = .0008067.*Sali
    TS = (0.14./96.062).*(Sal./1.80655); % in mol/kg-SW

    % CalculateK0:
    % Weiss, R. F., Marine Chemistry 2:203-215, 1974.
    TempK100  = TempK./100;
    lnK0 = -60.2409 + 93.4517 ./ TempK100 + 23.3585 .* log(TempK100) + Sal .*...
        (0.023517 - 0.023656 .* TempK100 + 0.0047036 .* TempK100 .^2);
    K0   = exp(lnK0);                  % this is in mol/kg-SW/atm

    % CalculateIonS:
    % This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4:
    IonS         = 19.924 .* Sal ./ (1000 - 1.005   .* Sal);

    % CalculateKS:

    switch WhoseKSO4
        case {1,3}
            % Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
            % The goodness of fit is .021.
            % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
            % TYPO on p. 121: the constant e9 should be e8.
            % This is from eqs 22 and 23 on p. 123, and Table 4 on p 121:
            lnKS = -4276.1./TempK + 141.328 - 23.093.*logTempK +...
            (-13856./TempK + 324.57 - 47.986.*logTempK).*sqrt(IonS) +...
                (35474./TempK - 771.54 + 114.723.*logTempK).*IonS +...
                (-2698./TempK).*sqrt(IonS).*IonS + (1776./TempK).*IonS.^2;

            KS = exp(lnKS)...            % this is on the free pH scale in mol/kg-H2O
                .* volkg;   % convert to mol/kg-SW
        case {2,4}
            % Khoo, et al, Analytical Chemistry, 49(1):29-34, 1977
            % KS was found by titrations with a hydrogen electrode
            % of artificial seawater containing sulfate (but without F)
            % at 3 salinities from 20 to 45 and artificial seawater NOT
            % containing sulfate (nor F) at 16 salinities from 15 to 45,
            % both at temperatures from 5 to 40 deg C.
            % KS is on the Free pH scale (inherently so).
            % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
            % He finds log(beta) which = my pKS;
            % his beta is an association constant.
            % The rms error is .0021 in pKS, or about .5% in KS.
            % This is equation 20 on p. 33:
            pKS = 647.59 ./ TempK - 6.3451 + 0.019085.*TempK - 0.5208.*sqrt(IonS);
            KS = 10.^(-pKS)...          % this is on the free pH scale in mol/kg-H2O
                .* volkg;    % convert to mol/kg-SW
        otherwise
            error('invalid WhoseKSO4 %g',WhoseKSO4);
    end

    % CalculateKF:
    % Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979:
    lnKF = 1590.2./TempK - 12.641 + 1.525.*IonS.^0.5;
    KF   = exp(lnKF)...                 % this is on the free pH scale in mol/kg-H2O
        .*volkg;          % convert to mol/kg-SW
    % Another expression exists for KF: Perez and Fraga 1987. Not used here since ill defined for low salinity. (to be used for S: 10-40, T: 9-33)
    % Nonetheless, P&F87 might actually be better than the fit of D&R79 above, which is based on only three salinities: [0 26.7 34.6]
    % lnKF = 874./TempK - 9.68 + 0.111.*Sal.^0.5; 
    % KF   = exp(lnKF);                   % this is on the free pH scale in mol/kg-SW

    % CalculatepHScaleConversionFactors:
    %       These are NOT pressure-corrected.
    SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);
    FREEtoTOT =  1 + TS./KS;

    % CalculatefH
    % Use GEOSECS's value for cases 1,2,3,4,5 (and 6) to convert pH scales.
    switch WhichKs
        case 8        
            fH = ones(nd); % this shouldn't occur in the program for this case
        case 7

            fH = 1.29 - 0.00204.*  TempK + (0.00046 -...
                0.00000148.*TempK).*Sal.*Sal;
            % Peng et al, Tellus 39B:439-458, 1987:
            % They reference the GEOSECS report, but round the value
            % given there off so that it is about .008 (1%) lower. It
            % doesn't agree with the check value they give on p. 456.
        otherwise

            fH = 1.2948 - 0.002036.*TempK + (0.0004607 -...
                0.000001475.*TempK).*Sal.^2;
            % Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition,
            % v. 3, 1982 (p. 80);
    end

    % CalculateKB:
    switch WhichKs
        case 8
            % Pure water case
            KB = zeros(nd);
        case {6,7}
            % This is for GEOSECS and Peng et al.
            % Lyman, John, UCLA Thesis, 1957
            % fit by Li et al, JGR 74:5507-5525, 1969:
            logKB = -9.26 + 0.00886.*Sal + 0.01.*TempC;
            KB = 10.^(logKB)...  % this is on the NBS scale
                ./fH;               % convert to the SWS scale
        otherwise
            % Dickson, A. G., Deep-Sea Research 37:755-766, 1990:
            lnKBtop = -8966.9 - 2890.53.*sqrSal - 77.942.*Sal +...
                1.728.*sqrSal.*Sal - 0.0996.*Sal.^2;
            lnKB = lnKBtop./TempK + 148.0248 + 137.1942.*sqrSal +...
                1.62142.*Sal + (-24.4344 - 25.085.*sqrSal - 0.2474.*...
                Sal).*logTempK + 0.053105.*sqrSal.*TempK;
            KB = exp(lnKB)...    % this is on the total pH scale in mol/kg-SW
                ./SWStoTOT;         % convert to SWS pH scale
    end

    % CalculateKW:

    switch WhichKs
        case 7
            % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
            lnKW = 148.9802 - 13847.26./TempK - 23.6521.*logTempK +...
                (-79.2447 + 3298.72./TempK + 12.0408.*logTempK).*...
                sqrSal - 0.019813.*Sal;
            KW = exp(lnKW); % this is on the SWS pH scale in (mol/kg-SW)^2
        case 8
            % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
            % refit data of Harned and Owen, The Physical Chemistry of
            % Electrolyte Solutions, 1958
            lnKW = 148.9802 - 13847.26./TempK - 23.6521.*logTempK;
            KW = exp(lnKW); % this is on the SWS pH scale in (mol/kg-SW)^2
        case 6
            KW = 0; % GEOSECS doesn't include OH effects
        otherwise
            % Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995.
            % his check value of 1.6 umol/kg-SW should be 6.2
            lnKW = 148.9802 - 13847.26./TempK - 23.6521.*logTempK +...
                (-5.977 + 118.67./TempK + 1.0495.*logTempK).*...
                sqrSal - 0.01615.*Sal;
            KW = exp(lnKW); % this is on the SWS pH scale in (mol/kg-SW)^2        
    end

    % CalculateKP1KP2KP3KSi:

    switch WhichKs
        case 7
            KP1 = 0.02*ones(nd);
            % Peng et al don't include the contribution from this term,
            % but it is so small it doesn't contribute. It needs to be
            % kept so that the routines work ok.
            % KP2, KP3 from Kester, D. R., and Pytkowicz, R. M.,
            % Limnology and Oceanography 12:243-252, 1967:
            % these are only for sals 33 to 36 and are on the NBS scale
            KP2 = exp(-9.039 - 1450./TempK)... % this is on the NBS scale
                ./fH;                          % convert to SWS scale
            KP3 = exp(4.466 - 7276./TempK)...  % this is on the NBS scale
                ./fH;                          % convert to SWS scale
            % Sillen, Martell, and Bjerrum,  Stability Constants of metal-ion complexes,
            % The Chemical Society (London), Special Publ. 17:751, 1964:
            KSi = 0.0000000004...              % this is on the NBS scale
                ./fH;                          % convert to SWS scale
        case {6,8}
            KP1 = zeros(nd); KP2 = zeros(nd); KP3 =zeros(nd); KSi = zeros(nd);
            % Neither the GEOSECS choice nor the freshwater choice
            % include contributions from phosphate or silicate.
        otherwise
            % Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
            % KP1, KP2, KP3 are on the SWS pH scale in mol/kg-SW.
            % KSi was given on the SWS pH scale in molal units.
            lnKP1 = -4576.752./TempK + 115.54 - 18.453.*logTempK + (-106.736./TempK +...
                0.69171).*sqrSal + (-0.65643./TempK - 0.01844).*Sal;
            KP1 = exp(lnKP1);
            lnKP2 = -8814.715./TempK + 172.1033 - 27.927.*logTempK + (-160.34./TempK +...
                1.3566).*sqrSal + (0.37335./TempK - 0.05778).*Sal;
            KP2 = exp(lnKP2);
            lnKP3 = -3070.75./TempK - 18.126 + (17.27039./TempK + 2.81197).*sqrSal +...
                (-44.99486./TempK - 0.09984).*Sal;
            KP3 = exp(lnKP3);
            lnKSi = -8904.2./TempK + 117.4 - 19.334.*logTempK + (-458.79./TempK +...
                3.5913).*sqrt(IonS) + (188.74./TempK - 1.5998).*IonS +...
                (-12.1652./TempK + 0.07871).*IonS.^2;
            KSi = exp(lnKSi)...                % this is on the SWS pH scale in mol/kg-H2O
                .*volkg;        % convert to mol/kg-SW
    end

    % CalculateK1K2:

    switch WhichKs
        case 1

            % ROY et al, Marine Chemistry, 44:249-267, 1993
            % (see also: Erratum, Marine Chemistry 45:337, 1994
            % and Erratum, Marine Chemistry 52:183, 1996)
            % Typo: in the abstract on p. 249: in the eq. for lnK1* the
            % last term should have S raised to the power 1.5.
            % They claim standard deviations (p. 254) of the fits as
            % .0048 for lnK1 (.5% in K1) and .007 in lnK2 (.7% in K2).
            % They also claim (p. 258) 2s precisions of .004 in pK1 and
            % .006 in pK2. These are consistent, but Andrew Dickson
            % (personal communication) obtained an rms deviation of about
            % .004 in pK1 and .003 in pK2. This would be a 2s precision
            % of about 2% in K1 and 1.5% in K2.
            % T:  0-45  S:  5-45. Total Scale. Artificial sewater.
            % This is eq. 29 on p. 254 and what they use in their abstract:
            lnK1 = 2.83655 - 2307.1266./TempK - 1.5529413.*logTempK +...
                (-0.20760841 - 4.0484./TempK).*sqrSal + 0.08468345.*Sal -...
                0.00654208.*sqrSal.*Sal;
            K1 = exp(lnK1)...               % this is on the total pH scale in mol/kg-H2O
                .*volkg...   % convert to mol/kg-SW
                ./SWStoTOT;                 % convert to SWS pH scale
            % This is eq. 30 on p. 254 and what they use in their abstract:
            lnK2 = -9.226508 - 3351.6106./TempK - 0.2005743.*logTempK +...
                (-0.106901773 - 23.9722./TempK).*sqrSal + 0.1130822.*Sal -...
                0.00846934.*sqrSal.*Sal;
            K2 = exp(lnK2)...               % this is on the total pH scale in mol/kg-H2O
                .*volkg...   % convert to mol/kg-SW
                ./SWStoTOT;                 % convert to SWS pH scale
        case 2
            % GOYET AND POISSON, Deep-Sea Research, 36(11):1635-1654, 1989
            % The 2s precision in pK1 is .011, or 2.5% in K1.
            % The 2s precision in pK2 is .02, or 4.5% in K2.
            % This is in Table 5 on p. 1652 and what they use in the abstract:
            pK1 = 812.27./TempK + 3.356 - 0.00171.*Sal.*logTempK...
                + 0.000091.*Sal.^2;
            K1 = 10.^(-pK1); % this is on the SWS pH scale in mol/kg-SW
            %
            % This is in Table 5 on p. 1652 and what they use in the abstract:
            pK2 = 1450.87./TempK + 4.604 - 0.00385.*Sal.*logTempK...
                + 0.000182.*Sal.^2;
            K2 = 10.^(-pK2); % this is on the SWS pH scale in mol/kg-SW
        case 3
            % HANSSON refit BY DICKSON AND MILLERO
            % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973
            % and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973.
            % on the SWS pH scale in mol/kg-SW.
            % Hansson gave his results on the Total scale (he called it
            % the seawater scale) and in mol/kg-SW.
            % Typo in DM on p. 1739 in Table 4: the equation for pK2*
            % for Hansson should have a .000132 *S^2
            % instead of a .000116 *S^2.
            % The 2s precision in pK1 is .013, or 3% in K1.
            % The 2s precision in pK2 is .017, or 4.1% in K2.
            % This is from Table 4 on p. 1739.
            pK1 = 851.4./TempK + 3.237 - 0.0106.*Sal + 0.000105.*Sal.^2;
            K1 = 10.^(-pK1); % this is on the SWS pH scale in mol/kg-SW
            %
            % This is from Table 4 on p. 1739.
            pK2 = -3885.4./TempK + 125.844 - 18.141.*logTempK...
                - 0.0192.*Sal + 0.000132.*Sal.^2;
            K2 = 10.^(-pK2); % this is on the SWS pH scale in mol/kg-SW
        case 4
            % MEHRBACH refit BY DICKSON AND MILLERO
            % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
            % on the SWS pH scale in mol/kg-SW.
            % Mehrbach et al gave results on the NBS scale.
            % The 2s precision in pK1 is .011, or 2.6% in K1.
            % The 2s precision in pK2 is .020, or 4.6% in K2.
            % Valid for salinity 20-40.
            % This is in Table 4 on p. 1739.
            pK1 = 3670.7./TempK - 62.008 + 9.7944.*logTempK...
                - 0.0118.*Sal + 0.000116.*Sal.^2;
            K1 = 10.^(-pK1); % this is on the SWS pH scale in mol/kg-SW
            %
            % This is in Table 4 on p. 1739.
            pK2 = 1394.7./TempK + 4.777 - 0.0184.*Sal + 0.000118.*Sal.^2;
            K2 = 10.^(-pK2); % this is on the SWS pH scale in mol/kg-SW
        case 5
            % HANSSON and MEHRBACH refit BY DICKSON AND MILLERO
            % Dickson and Millero, Deep-Sea Research,34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973,
            % Hansson, Acta Chemica Scandanavia, 27:931-944, 1973,
            % and Mehrbach et al, Limnol. Oceanogr.,18(6):897-907, 1973
            % on the SWS pH scale in mol/kg-SW.
            % Typo in DM on p. 1740 in Table 5: the second equation
            % should be pK2* =, not pK1* =.
            % The 2s precision in pK1 is .017, or 4% in K1.
            % The 2s precision in pK2 is .026, or 6% in K2.
            % Valid for salinity 20-40.
            % This is in Table 5 on p. 1740.
            pK1 = 845./TempK + 3.248 - 0.0098.*Sal + 0.000087.*Sal.^2;
            K1 = 10.^(-pK1); % this is on the SWS pH scale in mol/kg-SW
            %
            % This is in Table 5 on p. 1740.
            pK2 = 1377.3./TempK + 4.824 - 0.0185.*Sal + 0.000122.*Sal.^2;
            K2 = 10.^(-pK2); % this is on the SWS pH scale in mol/kg-SW
        case {6,7}
            % GEOSECS and Peng et al use K1, K2 from Mehrbach et al,
            % Limnology and Oceanography, 18(6):897-907, 1973.
            % I.e., these are the original Mehrbach dissociation constants.
            % The 2s precision in pK1 is .005, or 1.2% in K1.
            % The 2s precision in pK2 is .008, or 2% in K2.
            pK1 = - 13.7201 + 0.031334.*TempK + 3235.76./TempK...
                + 1.3e-5*Sal.*TempK - 0.1032.*Sal.^0.5;
            K1 = 10.^(-pK1)...         % this is on the NBS scale
                ./fH;                     % convert to SWS scale
            pK2 = 5371.9645 + 1.671221.*TempK + 0.22913.*Sal + 18.3802.*log10(Sal)...
                - 128375.28./TempK - 2194.3055.*log10(TempK) - 8.0944e-4.*Sal.*TempK...
                - 5617.11.*log10(Sal)./TempK + 2.136.*Sal./TempK; % pK2 is not defined for Sal=0, since log10(0)=-inf
            K2 = 10.^(-pK2)...         % this is on the NBS scale
                ./fH;                     % convert to SWS scale
        case 8
            % PURE WATER CASE
            % Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
            % K1 from refit data from Harned and Davis,
            % J American Chemical Society, 65:2030-2037, 1943.
            % K2 from refit data from Harned and Scholes,
            % J American Chemical Society, 43:1706-1709, 1941.
            % This is only to be used for Sal=0 water (note the absence of S in the below formulations)
            % These are the thermodynamic Constants:
            lnK1 = 290.9097 - 14554.21./TempK - 45.0575.*logTempK;
            K1 = exp(lnK1);
            lnK2 = 207.6548 - 11843.79./TempK - 33.6485.*logTempK;
            K2 = exp(lnK2);
        case 9
            % From Cai and Wang 1998, for estuarine use.
            % Data used in this work is from:
            % K1: Merhback (1973) for S>15, for S<15: Mook and Keone (1975)
            % K2: Merhback (1973) for S>20, for S<20: Edmond and Gieskes (1970)
            % Sigma of residuals between fits and above data: ±0.015, +0.040 for K1 and K2, respectively.
            % Sal 0-40, Temp 0.2-30
            % Limnol. Oceanogr. 43(4) (1998) 657-668
            % On the NBS scale
            % Their check values for F1 don't work out, not sure if this was correctly published...
            F1 = 200.1./TempK + 0.3220;
            pK1 = 3404.71./TempK + 0.032786.*TempK - 14.8435 - 0.071692.*F1.*Sal.^0.5 + 0.0021487.*Sal;
            K1  = 10.^-pK1...         % this is on the NBS scale
                ./fH;                    % convert to SWS scale (uncertain at low Sal due to junction potential);
            F2 = -129.24./TempK + 1.4381;
            pK2 = 2902.39./TempK + 0.02379.*TempK - 6.4980 - 0.3191.*F2.*Sal.^0.5 + 0.0198.*Sal;
            K2  = 10.^-pK2...         % this is on the NBS scale
                ./fH;                    % convert to SWS scale (uncertain at low Sal due to junction potential);
        case 10

            % From Lueker, Dickson, Keeling, 2000
            % This is Mehrbach's data refit after conversion to the total scale, for comparison with their equilibrator work.
            % Mar. Chem. 70 (2000) 105-119
            % Total scale and kg-sw
            pK1 = 3633.86./TempK-61.2172+9.6777.*log(TempK)-0.011555.*Sal+0.0001152.*Sal.^2;
            K1  = 10.^-pK1...           % this is on the total pH scale in mol/kg-SW
                ./SWStoTOT;                % convert to SWS pH scale
            pK2 = 471.78./TempK+25.929 -3.16967.*log(TempK)-0.01781 .*Sal+0.0001122.*Sal.^2;
            K2  = 10.^-pK2...           % this is on the total pH scale in mol/kg-SW
                ./SWStoTOT;                % convert to SWS pH scale
        case 11
            % Mojica Prieto and Millero 2002. Geochim. et Cosmochim. Acta. 66(14) 2529-2540.
            % sigma for pK1 is reported to be 0.0056
            % sigma for pK2 is reported to be 0.010
            % This is from the abstract and pages 2536-2537
            pK1 =  -43.6977 - 0.0129037.*Sal + 1.364e-4.*Sal.^2 + 2885.378./TempK +  7.045159.*log(TempK);
            pK2 = -452.0940 + 13.142162.*Sal - 8.101e-4.*Sal.^2 + 21263.61./TempK + 68.483143.*log(TempK)...
                + (-581.4428.*Sal + 0.259601.*Sal.^2)./TempK - 1.967035.*Sal.*log(TempK);
            K1 = 10.^-pK1; % this is on the SWS pH scale in mol/kg-SW
            K2 = 10.^-pK2; % this is on the SWS pH scale in mol/kg-SW
        case 12
            % Millero et al., 2002. Deep-Sea Res. I (49) 1705-1723.
            % Calculated from overdetermined WOCE-era field measurements
            % sigma for pK1 is reported to be 0.005
            % sigma for pK2 is reported to be 0.008
            % This is from page 1715
            pK1 =  6.359 - 0.00664.*Sal - 0.01322.*TempC + 4.989e-5.*TempC.^2;
            pK2 =  9.867 - 0.01314.*Sal - 0.01904.*TempC + 2.448e-5.*TempC.^2;
            K1 = 10.^-pK1; % this is on the SWS pH scale in mol/kg-SW
            K2 = 10.^-pK2; % this is on the SWS pH scale in mol/kg-SW
        case 13
            % From Millero 2006 work on pK1 and pK2 from titrations
            % Millero, Graham, Huang, Bustos-Serrano, Pierrot. Mar.Chem. 100 (2006) 80-94.
            % S=1 to 50, T=0 to 50. On seawater scale (SWS). From titrations in Gulf Stream seawater.
            pK1_0 = -126.34048 + 6320.813./TempK + 19.568224*log(TempK);
            A_1   = 13.4191*Sal.^0.5 + 0.0331.*Sal - 5.33e-5.*Sal.^2;
            B_1   = -530.123*Sal.^0.5 - 6.103.*Sal;
            C_1   = -2.06950.*Sal.^0.5;
            pK1= A_1 + B_1./TempK + C_1.*log(TempK) + pK1_0; % pK1 sigma = 0.0054
            K1 = 10.^-(pK1);
            pK2_0= -90.18333 + 5143.692./TempK + 14.613358*log(TempK);
            A_2   = 21.0894*Sal.^0.5 + 0.1248.*Sal - 3.687e-4.*Sal.^2;
            B_2   = -772.483*Sal.^0.5 - 20.051.*Sal;
            C_2   = -3.3336.*Sal.^0.5;
            pK2= A_2 + B_2./TempK + C_2.*log(TempK) + pK2_0; %pK2 sigma = 0.011
            K2 = 10.^-(pK2);
        case 14
            % From Millero, 2010, also for estuarine use.
            % Marine and Freshwater Research, v. 61, p. 139?142.
            % Fits through compilation of real seawater titration results:
            % Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), Millero et al. (2006)
            % Constants for K's on the SWS;
            % This is from page 141
            pK10 = -126.34048 + 6320.813./TempK + 19.568224.*log(TempK);
            % This is from their table 2, page 140.
            A1 = 13.4038.*Sal.^0.5 + 0.03206.*Sal - 5.242e-5.*Sal.^2;
            B1 = -530.659.*Sal.^0.5 - 5.8210.*Sal;
            C1 = -2.0664*Sal.^0.5;
            pK1 = pK10 + A1 + B1./TempK + C1.*log(TempK);
            K1 = 10.^-pK1;
            % This is from page 141
            pK20 =  -90.18333 + 5143.692./TempK + 14.613358.*log(TempK);
            % This is from their table 3, page 140.
            A2 = 21.3728.*Sal.^0.5 + 0.1218.*Sal - 3.688e-4.*Sal.^2;
            B2 = -788.289.*Sal.^0.5 - 19.189.*Sal;
            C2 = -3.374.*Sal.^0.5;
            pK2 = pK20 + A2 + B2./TempK + C2.*log(TempK);
            K2 = 10.^-pK2;
        otherwise
            error('unrecognized WhichKs %g',WhichKs);
    end

    % SD 2014-7-28 Calculate KH2S
    % Original data from
    % Millero, F. J., Plese, T., & Fernandez, M. (1988). Limnol. Oceanogr, 33(2), 269?274.
    % number here from
    % Hofmann AF etal  (2010)  Aquatic Geochemistry 16, 507?546, doi:10.1007/s10498-009-9084-1
    lnKH2S= (225.838 + 0.3449.*sqrSal - 0.0274*Sal) - 13275.3./TempK - 34.6435.*logTempK;
    KH2S = exp(lnKH2S)./SWStoTOT;         % convert to SW pH scale

    % SD 2015-01-03 Calculate KNH4
    % numbers from 
    % Hofmann AF etal  (2010)  Aquatic Geochemistry 16, 507?546, doi:10.1007/s10498-009-9084-1
    lnKNH4= (-0.25444 + 0.46532.*sqrSal - 0.01992*Sal) + (-6285.33-123.7184.*sqrSal+3.17556*Sal)./TempK + 0.0001635.*TempK;
    KNH4 = exp(lnKNH4);          % parameters are for SW pH scale

    %***************************************************************************
    %CorrectKsForPressureNow:
    % Currently: For WhichKs% = 1 to 7, all Ks (except KF and KS, which are on
    %       the free scale) are on the SWS scale.
    %       For WhichKs% = 6, KW set to 0, KP1, KP2, KP3, KSi don't matter.
    %       For WhichKs% = 8, K1, K2, and KW are on the "pH" pH scale
    %       (the pH scales are the same in this case); the other Ks don't matter.
    %
    %
    % No salinity dependence is given for the pressure coefficients here.
    % It is assumed that the salinity is at or very near Sali = 35.
    % These are valid for the SWS pH scale, but the difference between this and
    % the total only yields a difference of .004 pH units at 1000 bars, much
    % less than the uncertainties in the values.
    %****************************************************************************
    % The sources used are:
    % Millero, 1995:
    %       Millero, F. J., Thermodynamics of the carbon dioxide system in the
    %       oceans, Geochemica et Cosmochemica Acta 59:661-677, 1995.
    %       See table 9 and eqs. 90-92, p. 675.
    %       TYPO: a factor of 10^3 was left out of the definition of Kappa
    %       TYPO: the value of R given is incorrect with the wrong units
    %       TYPO: the values of the a's for H2S and H2O are from the 1983
    %                values for fresh water
    %       TYPO: the value of a1 for B(OH)3 should be +.1622
    %        Table 9 on p. 675 has no values for Si.
    %       There are a variety of other typos in Table 9 on p. 675.
    %       There are other typos in the paper, and most of the check values
    %       given don't check.
    % Millero, 1992:
    %       Millero, Frank J., and Sohn, Mary L., Chemical Oceanography,
    %       CRC Press, 1992. See chapter 6.
    %       TYPO: this chapter has numerous typos (eqs. 36, 52, 56, 65, 72,
    %               79, and 96 have typos).
    % Millero, 1983:
    %       Millero, Frank J., Influence of pressure on chemical processes in
    %       the sea. Chapter 43 in Chemical Oceanography, eds. Riley, J. P. and
    %       Chester, R., Academic Press, 1983.
    %       TYPO: p. 51, eq. 94: the value -26.69 should be -25.59
    %       TYPO: p. 51, eq. 95: the term .1700t should be .0800t
    %       these two are necessary to match the values given in Table 43.24
    % Millero, 1979:
    %       Millero, F. J., The thermodynamics of the carbon dioxide system
    %       in seawater, Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
    %       See table 5 and eqs. 7, 7a, 7b on pp. 1656-1657.
    % Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
    %       TYPO: the pressure dependence of K2 should have a 16.4, not 26.4
    %       This matches the GEOSECS results and is in Edmond and Gieskes.
    % Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic acid,
    %       boric acid, and the pH of seawater, Limnology and Oceanography
    %       13:403-417, 1968.
    % Edmond, John M. and Gieskes, J. M. T. M., The calculation of the degree of
    %       seawater with respect to calcium carbonate under in situ conditions,
    %       Geochemica et Cosmochemica Acta, 34:1261-1291, 1970.
    %****************************************************************************
    % These references often disagree and give different fits for the same thing.
    % They are not always just an update either; that is, Millero, 1995 may agree
    %       with Millero, 1979, but differ from Millero, 1983.
    % For WhichKs% = 7 (Peng choice) I used the same factors for KW, KP1, KP2,
    %       KP3, and KSi as for the other cases. Peng et al didn't consider the
    %       case of P different from 0. GEOSECS did consider pressure, but didn't
    %       include Phos, Si, or OH, so including the factors here won't matter.
    % For WhichKs% = 8 (freshwater) the values are from Millero, 1983 (for K1, K2,
    %       and KW). The other aren't used (TB = TS = TF = TP = TSi = 0.), so
    %       including the factors won't matter.
    %****************************************************************************
    %       deltaVs are in cm3/mole
    %       Kappas are in cm3/mole/bar
    %****************************************************************************

    %CorrectK1K2KBForPressure:

    switch WhichKs
        case 8

            %***PressureEffectsOnK1inFreshWater:
            %               This is from Millero, 1983.
            deltaV  = -30.54 + 0.1849 .*TempC - 0.0023366.*TempC.^2;
            Kappa   = (-6.22 + 0.1368 .*TempC - 0.001233 .*TempC.^2)./1000;
            lnK1fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
            %***PressureEffectsOnK2inFreshWater:
            %               This is from Millero, 1983.
            deltaV  = -29.81 + 0.115.*TempC - 0.001816.*TempC.^2;
            Kappa   = (-5.74 + 0.093.*TempC - 0.001896.*TempC.^2)./1000;
            lnK2fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
            lnKBfac = 0 ;%; this doesn't matter since TB = 0 for this case
        case {6,7}

            %               GEOSECS Pressure Effects On K1, K2, KB (on the NBS scale)
            %               Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982 quotes
            %               Culberson and Pytkowicz, L and O 13:403-417, 1968:
            %               but the fits are the same as those in
            %               Edmond and Gieskes, GCA, 34:1261-1291, 1970
            %               who in turn quote Li, personal communication
            lnK1fac = (24.2 - 0.085.*TempC).*Pbar./RT;
            lnK2fac = (16.4 - 0.04 .*TempC).*Pbar./RT;
            %               Takahashi et al had 26.4, but 16.4 is from Edmond and Gieskes
            %               and matches the GEOSECS results
            lnKBfac = (27.5 - 0.095.*TempC).*Pbar./RT;
        otherwise
            %***PressureEffectsOnK1:
            %               These are from Millero, 1995.
            %               They are the same as Millero, 1979 and Millero, 1992.
            %               They are from data of Culberson and Pytkowicz, 1968.
            deltaV  = -25.5 + 0.1271.*TempC;
            %                 'deltaV = deltaV - .151.*(Sali - 34.8); % Millero, 1979
            Kappa   = (-3.08 + 0.0877.*TempC)./1000;
            %                 'Kappa = Kappa  - .578.*(Sali - 34.8)/1000.; % Millero, 1979
            lnK1fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
            %               The fits given in Millero, 1983 are somewhat different.

            %***PressureEffectsOnK2:
            %               These are from Millero, 1995.
            %               They are the same as Millero, 1979 and Millero, 1992.
            %               They are from data of Culberson and Pytkowicz, 1968.
            deltaV  = -15.82 - 0.0219.*TempC;
            %                  'deltaV = deltaV + .321.*(Sali - 34.8); % Millero, 1979
            Kappa   = (1.13 - 0.1475.*TempC)./1000;
            %                 'Kappa = Kappa - .314.*(Sali - 34.8)./1000: % Millero, 1979
            lnK2fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
            %               The fit given in Millero, 1983 is different.
            %               Not by a lot for deltaV, but by much for Kappa. %

            %***PressureEffectsOnKB:
            %               This is from Millero, 1979.
            %               It is from data of Culberson and Pytkowicz, 1968.
            deltaV  = -29.48 + 0.1622.*TempC - 0.002608.*TempC.^2;
            %               Millero, 1983 has:
            %                 'deltaV = -28.56 + .1211.*TempCi - .000321.*TempCi.*TempCi
            %               Millero, 1992 has:
            %                 'deltaV = -29.48 + .1622.*TempCi + .295.*(Sali - 34.8)
            %               Millero, 1995 has:
            %                 'deltaV = -29.48 - .1622.*TempCi - .002608.*TempCi.*TempCi
            %                 'deltaV = deltaV + .295.*(Sali - 34.8); % Millero, 1979
            Kappa   = -2.84./1000; % Millero, 1979
            %               Millero, 1992 and Millero, 1995 also have this.
            %                 'Kappa = Kappa + .354.*(Sali - 34.8)./1000: % Millero,1979
            %               Millero, 1983 has:
            %                 'Kappa = (-3 + .0427.*TempCi)./1000
            lnKBfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
    end

    % CorrectKWForPressure:

    switch WhichKs
        case 8
            % PressureEffectsOnKWinFreshWater:
            %               This is from Millero, 1983.
            deltaV  =  -25.6 + 0.2324.*TempC - 0.0036246.*TempC.^2;
            Kappa   = (-7.33 + 0.1368.*TempC - 0.001233 .*TempC.^2)./1000;
            lnKWfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
            %               NOTE the temperature dependence of KappaK1 and KappaKW
            %               for fresh water in Millero, 1983 are the same.
        otherwise
            % GEOSECS doesn't include OH term, so this won't matter.
            % Peng et al didn't include pressure, but here I assume that the KW correction
            %       is the same as for the other seawater cases.
            % PressureEffectsOnKW:
            %               This is from Millero, 1983 and his programs CO2ROY(T).BAS.
            deltaV  = -20.02 + 0.1119.*TempC - 0.001409.*TempC.^2;
            %               Millero, 1992 and Millero, 1995 have:
            Kappa   = (-5.13 + 0.0794.*TempC)./1000; % Millero, 1983
            %               Millero, 1995 has this too, but Millero, 1992 is different.
            lnKWfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
            %               Millero, 1979 does not list values for these.
    end

    % PressureEffectsOnKF:
    %       This is from Millero, 1995, which is the same as Millero, 1983.
    %       It is assumed that KF is on the free pH scale.
    deltaV = -9.78 - 0.009.*TempC - 0.000942.*TempC.^2;
    Kappa = (-3.91 + 0.054.*TempC)./1000;
    lnKFfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
    % PressureEffectsOnKS:
    %       This is from Millero, 1995, which is the same as Millero, 1983.
    %       It is assumed that KS is on the free pH scale.
    deltaV = -18.03 + 0.0466.*TempC + 0.000316.*TempC.^2;
    Kappa = (-4.53 + 0.09.*TempC)./1000;
    lnKSfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

    % CorrectKP1KP2KP3KSiForPressure:
    % These corrections don't matter for the GEOSECS choice (WhichKs% = 6) and
    %       the freshwater choice (WhichKs% = 8). For the Peng choice I assume
    %       that they are the same as for the other choices (WhichKs% = 1 to 5).
    % The corrections for KP1, KP2, and KP3 are from Millero, 1995, which are the
    %       same as Millero, 1983.
    % PressureEffectsOnKP1:
    deltaV = -14.51 + 0.1211.*TempC - 0.000321.*TempC.^2;
    Kappa  = (-2.67 + 0.0427.*TempC)./1000;
    lnKP1fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
    % PressureEffectsOnKP2:
    deltaV = -23.12 + 0.1758.*TempC - 0.002647.*TempC.^2;
    Kappa  = (-5.15 + 0.09  .*TempC)./1000;
    lnKP2fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
    % PressureEffectsOnKP3:
    deltaV = -26.57 + 0.202 .*TempC - 0.003042.*TempC.^2;
    Kappa  = (-4.08 + 0.0714.*TempC)./1000;
    lnKP3fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
    % PressureEffectsOnKSi:
    %  The only mention of this is Millero, 1995 where it is stated that the
    %    values have been estimated from the values of boric acid. HOWEVER,
    %    there is no listing of the values in the table.
    %    I used the values for boric acid from above.
    deltaV = -29.48 + 0.1622.*TempC - 0.002608.*TempC.^2;
    Kappa  = -2.84./1000;
    lnKSifac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

    % SD 2014-7-28 Correct KH2S for pressure
    % From Hofmann AF etal  (2010)  Aquatic Geochemistry 16, 507?546, doi:10.1007/s10498-009-9084-1
    %         a0      a1               a2
    deltaV = -14.80 + 0.0020.*TempC - 0.4000e-3.*TempC.^2;
    %         b0      b1               (b2)
    Kappa  =  (2.89   + 0.0540.*TempC)./1000;
    lnKH2Sfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

    % SD 2015-01-03 Correct KNH4 for pressure
    % From Hofmann AF etal  (2010)  Aquatic Geochemistry 16, 507?546, doi:10.1007/s10498-009-9084-1
    %         a0      a1               a2
    deltaV = -26.43 + 0.0889.*TempC - 0.9050e-3.*TempC.^2;
    %         b0      b1               (b2)
    Kappa  =  (-5.03   + 0.0814.*TempC)./1000;
    lnKNH4fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

    % CorrectKsForPressureHere:
    K1fac  = exp(lnK1fac);  K1  = K1 .*K1fac;
    K2fac  = exp(lnK2fac);  K2  = K2 .*K2fac;
    KWfac  = exp(lnKWfac);  KW  = KW .*KWfac;
    KBfac  = exp(lnKBfac);  KB  = KB .*KBfac;
    KFfac  = exp(lnKFfac);  KF  = KF .*KFfac;
    KSfac  = exp(lnKSfac);  KS  = KS .*KSfac;
    KP1fac = exp(lnKP1fac); KP1 = KP1.*KP1fac;
    KP2fac = exp(lnKP2fac); KP2 = KP2.*KP2fac;
    KP3fac = exp(lnKP3fac); KP3 = KP3.*KP3fac;
    KSifac = exp(lnKSifac); KSi = KSi.*KSifac;
    KH2Sfac= exp(lnKH2Sfac);KH2S= KH2S.*KH2Sfac;
    KNH4fac= exp(lnKNH4fac);KNH4= KNH4.*KNH4fac;

    % CorrectpHScaleConversionsForPressure:
    % fH has been assumed to be independent of pressure.
    SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);
    FREEtoTOT =  1 + TS./KS;

    %  The values KS and KF are already pressure-corrected, so the pH scale
    %  conversions are now valid at pressure.

    % FindpHScaleConversionFactor:
    % this is the scale they will be put on

    switch pHScale
        case 1
            %Total
            pHfactor = SWStoTOT;
        case 2
            %SWS, they are all on this now
            pHfactor = ones(nd);
        case 3
            %pHfree
            pHfactor = SWStoTOT./FREEtoTOT;
        case 4 %pHNBS
            pHfactor = fH;
        otherwise
            error('invalid pHScale %g',pHScale);
    end
    % ConvertFromSWSpHScaleToChosenScale:
    K1  = K1.* pHfactor; K2  = K2.* pHfactor;
    KW  = KW.* pHfactor; KB  = KB.* pHfactor;
    KP1 = KP1.*pHfactor; KP2 = KP2.*pHfactor;
    KP3 = KP3.*pHfactor; KSi = KSi.*pHfactor;
    KH2S= KH2S.*pHfactor; KNH4 = KNH4.*pHfactor;

%     % CalculateFugacityConstants:
%     switch WhichKs
%         case {6,7}
%             % GEOSECS and Peng assume pCO2 = fCO2, or FugFac = 1
%             FugFac = ones(nd);
%         otherwise
%             % This assumes that the pressure is at one atmosphere, or close to it.
%             % Otherwise, the Pres term in the exponent affects the results.
%             %       Weiss, R. F., Marine Chemistry 2:203-215, 1974.
%             %       Delta and B in cm3/mol  
%             Delta = (57.7 - 0.118.*TempK);
%             b = -1636.75 + 12.0408.*TempK - 0.0327957.*TempK.^2 + 3.16528.*0.00001.*TempK.^3;
%             % For a mixture of CO2 and air at 1 atm (at low CO2 concentrations);
%             P1atm = 1.01325; % in bar
%             FugFac = exp((b + 2.*Delta).*P1atm./RT);
%     end

    % CalculateVPFac:
%     % Weiss, R. F., and Price, B. A., Nitrous oxide solubility in water and
%     %       seawater, Marine Chemistry 8:347-359, 1980.
%     % They fit the data of Goff and Gratch (1946) with the vapor pressure
%     %       lowering by sea salt as given by Robinson (1954).
%     % This fits the more complicated Goff and Gratch, and Robinson equations
%     %       from 273 to 313 deg K and 0 to 40 Sali with a standard error
%     %       of .015%, about 5 uatm over this range.
%     % This may be on IPTS-29 since they didn't mention the temperature scale,
%     %       and the data of Goff and Gratch came before IPTS-48.
%     % The references are:
%     % Goff, J. A. and Gratch, S., Low pressure properties of water from -160 deg
%     %       to 212 deg F, Transactions of the American Society of Heating and
%     %       Ventilating Engineers 52:95-122, 1946.
%     % Robinson, Journal of the Marine Biological Association of the U. K.
%     %       33:449-455, 1954.
%     %       This is eq. 10 on p. 350.
%     %       This is in atmospheres.
%     VPWP = exp(24.4543 - 67.4509.*(100./TempK) - 4.8489.*log(TempK./100));
%     VPCorrWP = exp(-0.000544.*Sal);
%     VPSWWP = VPWP.*VPCorrWP;
%     VPFac = 1 - VPSWWP; % this assumes 1 atmosphere

    % import results into C
%     TempK   = TempK;
%     logTempK = logTempK;
%     RT    = RT;
% 
%     C.TB  = TB;    % contemporary value from Sal
%     C.TF = TF;     % contemporary value from Sal
%     C.TS = TS;     % contemporary value
% 
%     C.fH    = fH;
%     C.FugFac = FugFac;
%     C.VPFac = VPFac;

    constants(:,I.K0)    = K0;   % carbonic acid
    constants(:,I.K1)    = K1;
    constants(:,I.K2)    = K2;
    constants(:,I.Kw)    = KW;   % water
    constants(:,I.Kb)    = KB;   % boric acid
    constants(:,I.KF)    = KF;   % hydrogen fluoride
    constants(:,I.Ks)    = KS;   % bisulfate
    constants(:,I.KP1)   = KP1;  % phosphoric acid
    constants(:,I.KP2)   = KP2;
    constants(:,I.KP3)   = KP3;
    constants(:,I.KSi)   = KSi;  % silicic acid
    
    switch WhichKs
    case {6,7}
        %
        % *** CalculateCaforGEOSECS:
        % Culkin, F, in Chemical Oceanography, ed. Riley and Skirrow, 1965:
        % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
        Ca = 0.01026.*Sal./35;
        % Culkin gives Ca = (.0213./40.078).*(Sal./1.80655) in mol/kg-SW
        % which corresponds to Ca = .01030.*Sal./35.
        %
        % *** CalculateKCaforGEOSECS:
        % Ingle et al, Marine Chemistry 1:295-307, 1973 is referenced in
        % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
        % but the fit is actually from Ingle, Marine Chemistry 3:301-319, 1975)
        KCa = 0.0000001.*(-34.452 - 39.866.*Sal.^(1./3) +...
            110.21.*log(Sal)./log(10) - 0.0000075752.*TempK.^2);
        % this is in (mol/kg-SW)^2
        %
        % *** CalculateKArforGEOSECS:
        % Berner, R. A., American Journal of Science 276:713-730, 1976:
        % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
        KAr = 1.45.*KCa;% ' this is in (mol/kg-SW)^2
        % Berner (p. 722) states that he uses 1.48.
        % It appears that 1.45 was used in the GEOSECS calculations
        %
        % *** CalculatePressureEffectsOnKCaKArGEOSECS:
        % Culberson and Pytkowicz, Limnology and Oceanography 13:403-417, 1968
        % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
        % but their paper is not even on this topic).
        % The fits appears to be new in the GEOSECS report.
        % I can't find them anywhere else.
        KCa = KCa.*exp((36   - 0.2 .*TempC).*Pbar./RT);
        KAr = KAr.*exp((33.3 - 0.22.*TempC).*Pbar./RT);
    otherwise
        % CalculateCa:
        % '       Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
        % '       this is .010285.*Sali./35
        Ca = 0.02128./40.087.*(Sal./1.80655);% ' in mol/kg-SW

        % CalciteSolubility:
        % '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
        logKCa = -171.9065 - 0.077993.*TempK + 2839.319./TempK;
        logKCa = logKCa + 71.595.*logTempK./log(10);
        logKCa = logKCa + (-0.77712 + 0.0028426.*TempK + 178.34./TempK).*sqrSal;
        logKCa = logKCa - 0.07711.*Sal + 0.0041249.*sqrSal.*Sal;
        % '       sd fit = .01 (for Sal part, not part independent of Sal)
        KCa = 10.^(logKCa);% ' this is in (mol/kg-SW)^2

        % AragoniteSolubility:
        % '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
        logKAr = -171.945 - 0.077993.*TempK + 2903.293./TempK;
        logKAr = logKAr + 71.595.*logTempK./log(10);
        logKAr = logKAr + (-0.068393 + 0.0017276.*TempK + 88.135./TempK).*sqrSal;
        logKAr = logKAr - 0.10018.*Sal + 0.0059415.*sqrSal.*Sal;
        % '       sd fit = .009 (for Sal part, not part independent of Sal)
        KAr    = 10.^(logKAr);% ' this is in (mol/kg-SW)^2

        % PressureCorrectionForCalcite:
        % '       Ingle, Marine Chemistry 3:301-319, 1975
        % '       same as in Millero, GCA 43:1651-1661, 1979, but Millero, GCA 1995
        % '       has typos (-.5304, -.3692, and 10^3 for Kappa factor)
        deltaVKCa = -48.76 + 0.5304.*TempC;
        KappaKCa  = (-11.76 + 0.3692.*TempC)./1000;
        lnKCafac  = (-deltaVKCa + 0.5.*KappaKCa.*Pbar).*Pbar./RT;
        KCa       = KCa.*exp(lnKCafac);

        % PressureCorrectionForAragonite:
        % '       Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979,
        % '       same as Millero, GCA 1995 except for typos (-.5304, -.3692,
        % '       and 10^3 for Kappa factor)
        deltaVKAr = deltaVKCa + 2.8;
        KappaKAr  = KappaKCa;
        lnKArfac  = (-deltaVKAr + 0.5.*KappaKAr.*Pbar).*Pbar./RT;
        KAr       = KAr.*exp(lnKArfac);
    end


    constants(:,I.Kcal)     =  KCa; % solubility products
    constants(:,I.Karg)     =  KAr;



end % end nested function

%%
function [constants] = calc_constants_OMG(I,T,D,S)

% calculate carbonate system equilibrium constants
    % using Dickson et al., (2007) Huide to Best Practices for Ocean CO2 Measurements
    % using Zeebe & Westbroeck (2001) CO2 in Seawater: Equilibrium, Kinectics, Isotopes
    % K1 & K2 are Merbach et al., 1973, refitted by Dickson & Millero 1987
    
    % set fit limits
    %T(T<2.0)=2.0;
    %T(T>35.0)=35.0;
    %S(S<26.0)=26.0;
    %S(S>43.0)=43.0;
    
    constants=zeros(size(T,1),13);
    
    % precalculate 
    TC = T;         % T in deg C
    T = T+273.15;   % T in Kelven
    Pbar=D./10;     % Pressure in bar (~depth)
    R=83.131;
    rRT=1./(R.*T);
    IonS = 19.924 .* S ./ (1000 - 1.005   .* S);
    conv_molarity_to_conc = (1 - 0.001005 .* S);
    
    % --------------------------------------------------------------------
    % Total F- & Total SO42- (mol/kg-SW)
    
    % from CO2SYS
    TF = (0.000067./18.998).*(S./1.80655);
    TS = (0.14./96.062).*(S./1.80655);
    
    % from GENIE
    %TF=0.00007.*S./35.0;
    %TS=0.02824.*S./35.0;
    
    % --------------------------------------------------------------------
    % KS (mol/kg-SW)
    
    KS = exp(...
         -4276.1./T + 141.328 - 23.093.*log(T) +...             
      (-13856./T + 324.57 - 47.986.*log(T)).*sqrt(IonS) +...     
      (35474./T - 771.54 + 114.723.*log(T)).*IonS +...           
      (-2698./T).*sqrt(IonS).*IonS + (1776./T).*IonS.^2);
  
    KS = KS .* conv_molarity_to_conc;
    
    % --------------------------------------------------------------------
    % KHF (mol/kg-SW)
    
    KF = exp( 1590.2./T - 12.641 + 1.525.*IonS.^0.5 );
    
    KF = KF .* conv_molarity_to_conc;
    
    % --------------------------------------------------------------------
    % Free pH to SWS conversion
    
    conv_SWStoTOT = (1 + TS./KS)./(1 + TS./KS + TF./KF);
    
    % --------------------------------------------------------------------
    % K1 (Merbach et al 1973, refit by Dickson and Millero 1987)
    constants(:,I.K1)=10.^-(3670.7./T-62.008+9.7944*log(T)...
        -0.0118.*S + 0.000116.*S.^2);
    
    delta_V=-25.50+0.1271.*TC+0.0e-3.*TC.^2;
    delta_ki=-3.08e-3+0.0877e-3.*TC;
    constants(:,I.K1)=exp((-delta_V+0.5.*delta_ki.*Pbar).*Pbar.*rRT).* constants(:,I.K1);
    
    % --------------------------------------------------------------------
	% K2 (Merbach et al 1973, refit by Dickson and Millero 1987)
    constants(:,I.K2)=10.^-(1394.7./T +4.777...
        - 0.0184.*S + 0.000118.*S.^2);
    
    delta_V=-15.82+-0.0219.*TC+0.0e-3.*TC.^2;
    delta_ki=-1.13e-3+-0.1475e-3.*TC;
    constants(:,I.K2)=exp((-delta_V+0.5.*delta_ki.*Pbar).*Pbar.*rRT).* constants(:,I.K2);

    % --------------------------------------------------------------------    
	% K0 (Weiss 1974)
	constants(:,I.K0)=exp( ...
	9345.17./T-60.2409+23.3585*log(T./100.0) ...
	+S.*(0.023517-0.00023656*T+0.0047036*(T./100.0).^2));
    
    % --------------------------------------------------------------------
	% KB (Dickson 1990b(
	constants(:,I.Kb)=exp( ...
	(-8966.90 - 2890.53.*S.^0.5 - 77.942.*S + 1.728.*S.^(3.0/2.0) - 0.0996.*S.^2.0)./T ...
	+ 148.0248 + 137.1942.*S.^0.5 + 1.62142.*S ...
	-(24.4344 + 25.085.*S.^0.5 + 0.2474.*S).*log(T) ...
	+ 0.053105.*S.^0.5.*T);

    constants(:,I.Kb) = constants(:,I.Kb) ./ conv_SWStoTOT;

    delta_V=-29.48 + 0.1622.*TC - 0.002608.*TC.^2;
    delta_ki=-2.84e-3 + 0.0e-3.*TC;
    constants(:,I.Kb)=exp((-delta_V+0.5.*delta_ki.*Pbar).*Pbar.*rRT).* constants(:,I.Kb);
        
    % --------------------------------------------------------------------
	% Kw (Millero 1995)
	constants(:,I.Kw)=exp( ...
	148.9802 - 13847.26./T  -23.6521.*log(T) ...
	+(118.67./T  -5.977 + 1.0495.*log(T)).*S.^0.5 - 0.01615.*S);

    delta_V= -20.02 + 0.1119.*TC + -0.001409.*TC.^2;
    delta_ki=-5.13e-3 + 0.0794e-3.*TC;
    constants(:,I.Kw)=exp((-delta_V+0.5.*delta_ki.*Pbar).*Pbar.*rRT).* constants(:,I.Kw);
    
    % --------------------------------------------------------------------
	%KSi (Millero 1995)
	loc_I=(19.924.*S)./(1000.0-1.005.*S);
	%I=0.02*S

	constants(:,I.KSi)=exp( ...
	-8904.2./T  +117.4 - 19.334.*log(T) ...
	+((-458.79./T + 3.5913).*(loc_I.^(-0.5)) ...
	+(188.74./T - 1.5998)).*loc_I ...
	+(-12.1652./T  +0.07871).*(loc_I.^2) ...
	+log(1.0 - 0.001005.*S));

    delta_V=-29.48 + 0.1622.*TC - 0.002608.*TC.^2;
    delta_ki=-0.0028+0.0e-3.*TC;
    constants(:,I.KSi)=exp((-delta_V+0.5.*delta_ki.*Pbar).*Pbar.*rRT).* constants(:,I.KSi);
    
    % --------------------------------------------------------------------
	% KP1 (Millero 1995)
	constants(:,I.KP1)=exp( ...
	-4576.752./T + 115.54 - 18.453.*log(T) ...
	+(-106.736./T + 0.69171).*S.^0.5 ...
	+(-0.65643./T - 0.01844).*S);

    delta_V=-14.51+0.1211.*TC+-0.321e-3.*TC.^2;
    delta_ki=-2.67e-3+0.0427e-3.*TC;
    constants(:,I.KP1)=exp((-delta_V+0.5.*delta_ki.*Pbar).*Pbar.*rRT).* constants(:,I.KP1);

    % --------------------------------------------------------------------
	% KP2 (Millero 1995)
	constants(:,I.KP2)=exp( ...
	-8814.715./T + 172.1033 - 27.927.*log(T) ...
	+(-160.34./T + 1.3566).*S.^0.5 ...
	+(0.37335./T  -0.05778).*S);

    delta_V=-23.12+0.1758.*TC+-2.647e-3.*TC.^2;
    delta_ki=-5.15e-3+0.09e-3.*TC;
    constants(:,I.KP2)=exp((-delta_V+0.5.*delta_ki.*Pbar).*Pbar.*rRT).* constants(:,I.KP2);

    % --------------------------------------------------------------------
	% KP3 (millero 1995)
	constants(:,I.KP3)=exp( ...
	-3070.75./T - 18.126 ...
	+(17.27039./T + 2.81197).*S.^0.5 ...
	+(-44.99486./T - 0.09984).*S);

    delta_V=-26.57+0.2020.*TC+-3.042e-3.*TC.^2;
    delta_ki=-4.08e-3+0.0714e-3.*TC;
    constants(:,I.KP3)=exp((-delta_V+0.5.*delta_ki.*Pbar).*Pbar.*rRT).* constants(:,I.KP3);

    % --------------------------------------------------------------------
    % Solubility product of Calcite (Mucci 1983)
    constants(:,I.Kcal)=10.^( ...
       -171.9065 - 0.077993.*T + 2839.319./T ...
       + 71.595 .* log10(T) ...
       +(-0.77712 + 0.0028426 .*T + 178.34 ./T) .* S.^0.5 ...
       - 0.07711 .*S + 0.0041249 .* S.^1.5);
   
    delta_V=-48.76+0.5304.*TC+0.0e-3.*TC.^2;
    delta_ki=-11.76e-3+0.3692e-3.*TC;
    constants(:,I.Kcal)=exp((-delta_V+0.5.*delta_ki.*Pbar).*Pbar.*rRT).* constants(:,I.Kcal);

   % --------------------------------------------------------------------    
   % Solubility product of Aragonite (Mucci 1983)
   constants(:,I.Karg)=10.^( ...
       -171.945 - 0.077993.*T + 2903.293./T ...
       + 71.595 .* log10(T) ...
       +(-0.068393 + 0.0017276 .*T + 88.135 ./T) .* S.^0.5 ...
       - 0.10018 .*S + 0.0059415 .* S.^1.5);
   
    delta_V=-45.96+0.5304.*TC+0.0e-3.*TC.^2;
    delta_ki=-11.76e-3+0.3692e-3.*TC;
    constants(:,I.Karg)=exp((-delta_V+0.5.*delta_ki.*Pbar).*Pbar.*rRT).* constants(:,I.Karg);
    
   % --------------------------------------------------------------------    
   % KS 
   constants(:,I.Ks)=KS;
   
   delta_V = -18.03 + 0.0466.*TC + 0.000316.*TC.^2;
   delta_ki = (-4.53 + 0.09.*TC)./1000;
   constants(:,I.Ks)=exp((-delta_V+0.5.*delta_ki.*Pbar).*Pbar.*rRT).* constants(:,I.Ks);
   
   % --------------------------------------------------------------------    
   % KF
   constants(:,I.KF)=KF;
   
   delta_V = -9.78 - 0.009.*TC - 0.000942.*TC.^2;
   delta_ki = (-3.91 + 0.054.*TC)./1000;
   constants(:,I.KF)=exp((-delta_V+0.5.*delta_ki.*Pbar).*Pbar.*rRT).* constants(:,I.KF);

end

%%
function [constants_interp] = interpolate_constants (OCEAN, I , ocn_pars , gchem_pars)

    constants_interp=zeros(ocn_pars.nb,13,ocn_pars.n_A);
    for n=1:ocn_pars.n_A
        constants_interp(:,:,n)=calc_constants(I,gchem_pars,OCEAN.T(:,n),OCEAN.S(:,n),ocn_pars.depth);
    end

end

%%
function [constants] = calc_constants(I,gchem_pars,T,S,D)

    if strcmp(gchem_pars.constants_method,'OMG')
        constants=calc_constants_OMG(I,T,D,S);
    elseif strcmp(gchem_pars.constants_method,'CO2SYS')
        constants=calc_constants_CO2SYS(I,T,D,S,gchem_pars.pHScale,gchem_pars,K1K2,gchem_pars.KSO4);
    end

end

%%
function [ECC] = calc_carbonate_constants(ECC,parameters,functions,full_calc,tcyc)

        if parameters.bgc_pars.CARBCHEM_select

        % initialise variables
        if ~full_calc
            ind=parameters.ocn_pars.Ib;
        else
            ind=[1:1:parameters.ocn_pars.nb]';
        end


        % calculate constants
        if parameters.gchem_pars.interpolate_constants
            for n=1:13
                ECC.constants(ind,n)=functions.OMG_fcns.interpolate_var_at_t(tcyc,parameters.ocn_pars,ECC.constants_interp(ind,n,:));
            end
        else
            ECC.constants(ind,:)=functions.gchem_fcns.calc_constants(parameters.ind_pars,parameters.gchem_pars,parameters.ocn_pars.T(ind),parameters.ocn_pars.S(ind),parameters.ocn_pars.depth(ind));
        end

        end

end

%%
function [ ECC ] = solve_carbonate_system (ECC,TRACERS,parameters,full_calc)

    % Follows et al., (2006) Ocean Modelling

    % dic = dissolved inorganic carbon; pt = dissolved inorganic phosphorus
    % sit = dissolved inorganic silica, bt = dissolved inorganic boron
    % ta = total alkalinity; ca = carbonate alkalinity; H = [H+]
    % pCO2 = partial pressure CO2; ff = fugacity of CO2
    % k1, k2 = carbonate equilibrium coeffs; kw = dissociation of water
    % klp, k2p, k3p = phosphate equilibrium coefficients
    % ksi, kb = silicate and borate equilibrium coefficients
    % Equilibrium relationships from DOE handbook (DOE, 1994):
    % coefficients evaluated elsewhere and passed in.
    
    if parameters.bgc_pars.CARBCHEM_select

        % initialise variables
        if ~full_calc
            ind=parameters.ocn_pars.Ib;
        else
            ind=[1:1:parameters.ocn_pars.nb]';
        end

        I=parameters.ind_pars;
        dt_ind=parameters.ocn_pars.Idt;

        sit=zeros(numel(TRACERS(ind,I.DIC)),1); % treated as per cGENIE
        bt=ones(numel(TRACERS(ind,I.DIC)),1).*(0.0004106*(parameters.ocn_pars.S(ind,1)./35.0)); % total boron conc. from ZW2001 (mol kg-1)
        % from CO2SYS
        %TF = (0.000067./18.998).*(s./1.80655);
        %TS = (0.14./96.062).*(s./1.80655);
        % from GENIE
        TF=0.00007.*parameters.ocn_pars.S(ind,1)./35.0;
        TS=0.02824.*parameters.ocn_pars.S(ind,1)./35.0;

        %Ca=0.01028*s./35.0;
        Ca=0.01028*parameters.ocn_pars.S(ind,1)./35.0;

        % First guess of [H+]: from last timestep
        rel_tol=1e-5; % relative tolerance for convergence of [H+]
        hg=ECC.state(ind,I.H); % incoming [H+]
        H=ones(numel(hg),1); % new guess (set to 1.0 initially)

        H_hg_rel_diff=abs(H-hg)./hg;
        n_loop=0;

        while any(H_hg_rel_diff>rel_tol)

            % estimate contributions to total alk from borate, silicate, phosphate
            bohg = bt .* ...
                ECC.constants(ind,I.Kb) ./ ...
                ( hg + ECC.constants(ind,I.Kb) );

            %siooh3g = sit.*ksi./(ksi + hg);
            siooh3g = sit .* ...
                ECC.constants(ind,I.KSi) ./ ...
                ( hg + ECC.constants(ind,I.KSi) );

            %denom = hg.*hg.*hg + (k1p.*hg.*hg) + (k1p.*k2p.*hg) + (k1p.*k2p.*k3p);
            denom = hg .* hg .* hg + ...
                ( ECC.constants(ind,I.KP1) .*hg .*hg ) + ...
                ( ECC.constants(ind,I.KP1) .* ECC.constants(ind,I.KP2) .* hg ) + ...
                ( ECC.constants(ind,I.KP1) .* ECC.constants(ind,I.KP2) .* ECC.constants(ind,I.KP3)) ;

            %h3po4g = (pt.*hg.*hg.*hg)./denom;
            h3po4g = ( TRACERS(ind,I.PO4) .* hg .* hg .* hg ) ./ denom;

            %h2po4g = (pt.*k1p.*hg.*hg)./denom;
            h2po4g = ( TRACERS(ind,I.PO4) .* ...
                ECC.constants(ind,I.KP1) .*hg .* hg ) ./ denom;

            %hpo4g = (pt.*k1p.*k2p.*hg)./denom;
            hpo4g = ( TRACERS(ind,I.PO4) .* ...
                ECC.constants(ind,I.KP1) .* ECC.constants(ind,I.KP2) .* hg ) ./ denom;

            %po4g = (pt.*k1p.*k2p.*k3p)./denom;
            po4g = ( TRACERS(ind,I.PO4) .* ...
                ECC.constants(ind,I.KP1) .* ECC.constants(ind,I.KP2) .* ECC.constants(ind,I.KP3) ) ./ denom;

            %hso4 = TS ./ (1.0 + ks./hg);

            %hf = TF ./ (1.0 + kf./hg);

            % estimate carbonate alkalinity
            fg = - bohg - (ECC.constants(ind,I.Kw)./hg) + hg - hpo4g - 2.0*po4g + h3po4g - siooh3g ;
            %fg = - bohg - (kw./hg) + hg - hpo4g - 2.0*po4g + h3po4g - siooh3g + hf - hso4;

            %cag = ta + fg;
            cag = TRACERS(ind,I.ALK) + fg;

            % improved estimate of hydrogen ion conc
            %gamm = dic./cag;
            gamm = TRACERS(ind,I.DIC) ./ cag;


            %dummy = (1.0-gamm).*(1.0-gamm).*k1.*k1 - 4.0.*k1.*k2.*(1.0 - 2.0.*gamm);
            dummy = (1.0 - gamm) .* (1.0-gamm) .* ...
                ECC.constants(ind,I.K1) .* ECC.constants(ind,I.K1) - ...
                4.0 .* ECC.constants(ind,I.K1) .* ECC.constants(ind,I.K2) .* ...
                (1.0 - 2.0.*gamm);

            %H = 0.5*((gamm-1.0).*k1 + sqrt(dummy));
            H = 0.5 * ((gamm-1.0) .* ECC.constants(ind,I.K1) + sqrt(dummy));

            H_hg_rel_diff=abs(H-hg)./hg; % relative difference between initial H+ and new guess
            hg=H; % update [H+]
            n_loop=n_loop+1;

            % stop OMG if [H+] does not converge to tolerance limit
            if n_loop>20
                error('Excessive iterations of carbonate chemistry')
                return
            end

        end

        % output variables
        % n.b. only certain ones here to save time!
        denom_2= H .* H + ECC.constants(ind,I.K1) .* H + ECC.constants(ind,I.K1) .* ECC.constants(ind,I.K2);

        ECC.state(ind,I.CO2) = TRACERS(ind,I.DIC) ./ (1.0 + (ECC.constants(ind,I.K1) ./ H ) + ...
            ( ECC.constants(ind,I.K1) .* ECC.constants(ind,I.K2) ./ (H .* H) ));

        ECC.state(ind,I.CO3) = TRACERS(ind,I.DIC) .* ECC.constants(ind,I.K1) .* ECC.constants(ind,I.K2) ./ denom_2;

        ECC.state(ind,I.SAT_CA) = Ca .* ECC.state(ind,I.CO3) ./ ECC.constants(ind,I.Kcal);

        ECC.state(ind,I.H)=H;

        % output everything if specified
        if full_calc
            ECC.state(ind,I.HCO3  ) = TRACERS(ind,I.DIC) .* ECC.constants(ind,I.K1) .* H ./ denom_2;
            ECC.state(ind,I.PH    ) = -log10(H);
            ECC.state(ind,I.SAT_AR) = Ca.*ECC.state(ind,I.CO3)./ECC.constants(ind,I.Karg);
        end
    else % (if parameters.bgc_pars.CARBCHEM_select not true)
        % pass out empty ECC.State array;
        ECC.state = [];
    end
    

end
