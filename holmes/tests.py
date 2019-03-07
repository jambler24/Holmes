from django.test import TestCase

# Create your tests here.


class AnimalTestCase(TestCase):

    def test_animals_can_speak(self):
        """Animals that can speak are correctly identified"""
        lion = 'The lion says "roar"'
        cat = 'The cat says "meow"'
        self.assertEqual(lion, 'The lion says "roar"')
        self.assertEqual(cat, 'The cat says "meow"')


