# convenient makefile

.PHONY: clean
clean:
	$(RM) *.pyc *.swp
	$(RM) -r  __pycache__/ .cache/
