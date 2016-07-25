from ..classes import Population, Cage, GATools
import pytest

def test_init():
    Population(Cage.__new__(Cage), Cage.__new__(Cage),
               Population.__new__(Population), Cage.__new__(Cage), 
               Population.__new__(Population), GATools.__new__(GATools))

    with pytest.raises(TypeError):
       Population(Cage.__new__(Cage), Cage.__new__(Cage), 
                  Population.__new__(Population), Cage.__new__(Cage),
                  Population.__new__(Population), 
                  GATools.__new__(GATools), GATools.__new__(GATools))

    with pytest.raises(TypeError):        
        Population(Cage.__new__(Cage), Cage.__new__(Cage), 
                   Population.__new__(Population), Cage.__new__(Cage),
                   Population.__new__(Population), 
                   GATools.__new__(GATools), 12)

    Population(*[None for x in range(0,10)])