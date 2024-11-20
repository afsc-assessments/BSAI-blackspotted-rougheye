if NOT exist .\results mkdir .\results   

copy rougheye24_%1.dat rougheye24.dat
rougheye24 -nox 

move rougheye24.rep .\results\m_ai_24_1_%1.rep
move rougheye24.std .\results\m_ai_24_1_%1.std
move rougheye24.par .\results\m_ai_24_1_%1.par
move rougheye24.rdat .\results\m_ai_24_1_%1.rdat
move rougheye24.cor  .\results\m_ai_24_1_%1.cor
